import os
import shutil
import numpy as np
import wandb
import torch
import torch.distributed as dist
import pytorch_lightning as pl
from pytorch_lightning import seed_everything
from pytorch_lightning.strategies import DDPStrategy
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataLoader
from torch_geometric.transforms import Compose
from tqdm.auto import tqdm

import hydra
from omegaconf import DictConfig

import utils.misc as misc
import utils.train as utils_train
import utils.transforms as trans
from datasets import get_dataset
from datasets.pl_data import FOLLOW_BATCH
from models.molopt_score_model import ScorePosNet3D
from models.pl_model import MoleculeTrainer

from graphbap.bapnet import BAPNet

import warnings
warnings.filterwarnings("ignore")


def load_pretrained_weights(model, ckpt_path):
    if ckpt_path in (None, '', 'null'):
        return
    ckpt_path = os.path.expanduser(str(ckpt_path))
    if not os.path.exists(ckpt_path):
        raise FileNotFoundError(f'pretrained checkpoint not found: {ckpt_path}')
    checkpoint = torch.load(ckpt_path, map_location='cpu', weights_only=False)
    state_dict = checkpoint.get('state_dict', checkpoint)
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    print(f'Loaded pretrained checkpoint: {ckpt_path}')
    print(f'Missing keys: {len(missing)}; unexpected keys: {len(unexpected)}')
    if missing:
        print('Missing key examples:', missing[:20])
    if unexpected:
        print('Unexpected key examples:', unexpected[:20])


@hydra.main(config_path='./configs', config_name='config')
def main(cfg: DictConfig):
    seed_everything(cfg.train.seed, workers=True)

    wandb_logger = misc.init_wandb(cfg)
    ckpt_dir, vis_dir = misc.create_folders(cfg)

    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_featurizer = trans.FeaturizeLigandAtom(cfg.data.transform.ligand_atom_mode)
    transform_list = [
        protein_featurizer,
        ligand_featurizer,
        trans.FeaturizeLigandBond(),
    ]
    if cfg.data.transform.random_rot:
        transform_list.append(trans.RandomRotation())
    transform = Compose(transform_list)

    subsets = get_dataset(
        config=cfg.data,
        transform=transform,
    )

    train_set, val_set, test_set = subsets['train'], subsets['valid'], subsets['test']
    print(f'Training: {len(train_set)} Validation: {len(val_set)} Testing: {len(test_set)}')

    collate_exclude_keys = ['ligand_nbh_list']
    num_workers = int(getattr(cfg.train, 'num_workers', 0))
    val_num_workers = int(getattr(cfg.train, 'val_num_workers', num_workers))
    pin_memory = bool(getattr(cfg.train, 'pin_memory', True))
    persistent_workers = bool(getattr(cfg.train, 'persistent_workers', True))
    prefetch_factor = int(getattr(cfg.train, 'prefetch_factor', 2))

    train_loader_kwargs = {
        'num_workers': num_workers,
        'pin_memory': pin_memory,
        'persistent_workers': persistent_workers and num_workers > 0,
    }
    val_loader_kwargs = {
        'num_workers': val_num_workers,
        'pin_memory': pin_memory,
        'persistent_workers': persistent_workers and val_num_workers > 0,
    }
    if num_workers > 0:
        train_loader_kwargs['prefetch_factor'] = prefetch_factor
    if val_num_workers > 0:
        val_loader_kwargs['prefetch_factor'] = prefetch_factor

    train_loader = DataLoader(
        train_set,
        batch_size=cfg.train.batch_size,
        shuffle=True,
        follow_batch=FOLLOW_BATCH,
        exclude_keys=collate_exclude_keys,
        **train_loader_kwargs,
    )
    val_loader = DataLoader(
        val_set,
        batch_size=cfg.train.val_batch_size,
        shuffle=False,
        follow_batch=FOLLOW_BATCH,
        exclude_keys=collate_exclude_keys,
        **val_loader_kwargs,
    )

    net_cond = BAPNet(ckpt_path=cfg.net_cond.ckpt_path,
                      hidden_nf=cfg.net_cond.hidden_dim)
    net_cond.freeze_the_model()

    pretrained_ckpt = getattr(cfg.train, 'pretrained_ckpt', None)
    if bool(getattr(cfg.train, 'freeze_backbone', False)) and pretrained_ckpt in (None, '', 'null'):
        raise RuntimeError(
            'train.freeze_backbone=true requires train.pretrained_ckpt to be set; '
            'otherwise you would freeze a randomly-initialised backbone and only '
            'train the gate on top of noise. Provide the Apo2Mol checkpoint path.'
        )

    model = MoleculeTrainer(
        cfg,
        protein_featurizer,
        ligand_featurizer,
        net_cond,
        )
    load_pretrained_weights(model, pretrained_ckpt)

    # Define the checkpoint callback
    checkpoint_callback = ModelCheckpoint(
        dirpath=ckpt_dir,
        filename='{epoch}',
        save_top_k=5,
        monitor='validation/loss',
        mode='min',
        save_last=True,
        every_n_epochs=int(getattr(cfg.train, 'val_freq', 1)),
    )

    early_stop_callback = EarlyStopping(
        monitor='validation/loss',
        patience=cfg.train.patience,
        verbose=True,
        mode='min'
    )

    # Use Tensor Cores
    torch.set_float32_matmul_precision('high')

    strategy = cfg.sys.strategy
    if cfg.sys.strategy == 'ddp':
        strategy = DDPStrategy(
            find_unused_parameters=bool(getattr(cfg.sys, 'find_unused_parameters', False)),
            gradient_as_bucket_view=bool(getattr(cfg.sys, 'gradient_as_bucket_view', True)),
        )

    trainer = pl.Trainer(
        accelerator=cfg.sys.accelerator,
        strategy=strategy,
        devices=cfg.sys.devices,
        num_nodes=cfg.sys.num_nodes,
        max_epochs=cfg.train.max_iters,
        max_steps=int(getattr(cfg.train, 'max_steps', -1)),
        limit_train_batches=getattr(cfg.train, 'limit_train_batches', 1.0),
        limit_val_batches=getattr(cfg.train, 'limit_val_batches', 1.0),
        check_val_every_n_epoch=int(getattr(cfg.train, 'val_freq', 1)),
        accumulate_grad_batches=int(getattr(cfg.train, 'n_acc_batch', 1)),
        num_sanity_val_steps=0,
        gradient_clip_val=cfg.train.max_grad_norm,
        precision=getattr(cfg.sys, 'precision', '32-true'),
        sync_batchnorm=bool(getattr(cfg.sys, 'sync_batchnorm', False)),
        benchmark=bool(getattr(cfg.sys, 'benchmark', True)),
        log_every_n_steps=int(getattr(cfg.train, 'log_every_n_steps', 50)),
        logger=wandb_logger,
        callbacks=[checkpoint_callback, early_stop_callback],
    )

    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)

if __name__ == '__main__':
    main()
