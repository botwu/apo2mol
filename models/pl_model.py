import numpy as np
import torch
import torch.distributed as dist
import pytorch_lightning as pl
from sklearn.metrics import roc_auc_score
from tqdm.auto import tqdm

import utils.train as utils_train
import utils.transforms as trans
from utils.data import apply_transforms_tensor_batch
from datasets.pl_data import FOLLOW_BATCH
from models.molopt_score_model import ScorePosNet3D

from graphbap.bapnet import BAPNet


def get_auroc(y_true, y_pred, feat_mode):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    avg_auroc = 0.
    total_weight = 0
    possible_classes = set(y_true)
    for c in possible_classes:
        binary_target = y_true == c
        if np.unique(binary_target).size < 2:
            continue
        auroc = roc_auc_score(binary_target, y_pred[:, c])
        class_weight = np.sum(binary_target)
        avg_auroc += auroc * class_weight
        total_weight += class_weight
        mapping = {
            'basic': trans.MAP_INDEX_TO_ATOM_TYPE_ONLY,
            'add_aromatic': trans.MAP_INDEX_TO_ATOM_TYPE_AROMATIC,
            'full': trans.MAP_INDEX_TO_ATOM_TYPE_FULL
        }
        print(f'atom: {mapping[feat_mode][c]} \t auc roc: {auroc:.4f}')
    return avg_auroc / total_weight if total_weight > 0 else float('nan')


class MoleculeTrainer(pl.LightningModule):
    def __init__(self, config, protein_featurizer, ligand_featurizer, net_cond):
        super(MoleculeTrainer, self).__init__()
        self.config = config
        self.protein_featurizer = protein_featurizer
        self.ligand_featurizer = ligand_featurizer
        self.net_cond = net_cond

        self.model = ScorePosNet3D(
            config.model,
            protein_atom_feature_dim=protein_featurizer.feature_dim,
            ligand_atom_feature_dim=ligand_featurizer.feature_dim
        )
        self._apply_freeze_policy()
        self.save_hyperparameters()

        # For training logging
        self.train_iterations = 0
        # For validation logging
        self.all_pred_v, self.all_true_v = [], []
        self.sum_loss, self.sum_loss_ligand_pos, self.sum_loss_v, self.sum_n = 0, 0, 0, 0
        self.sum_loss_protein_tr, self.sum_loss_protein_rot, self.sum_loss_protein_chi = 0, 0, 0
        self.sum_rmsd_protein_pos = 0

    def forward(self, *args, **kwargs):
        return self.model(*args, **kwargs)

    def _apply_freeze_policy(self):
        train_cfg = self.config.train
        freeze_backbone = bool(getattr(train_cfg, 'freeze_backbone', False))
        freeze_residue_head = bool(getattr(train_cfg, 'freeze_residue_head', True))
        if not freeze_backbone:
            return
        for name, param in self.model.named_parameters():
            if name.startswith('cross_attn_gate.'):
                param.requires_grad = True
                continue
            if not freeze_residue_head and name.startswith('res_inference.'):
                param.requires_grad = True
                continue
            param.requires_grad = False

    def configure_optimizers(self):
        opt_cfg = self.config.train.optimizer
        trainable = [p for p in self.model.parameters() if p.requires_grad]
        if len(trainable) == 0:
            raise RuntimeError(
                'No trainable parameters found. Check train.freeze_backbone and model.pocket_router_mode.'
            )
        if opt_cfg.type == 'adam':
            optimizer = torch.optim.Adam(
                trainable,
                lr=opt_cfg.lr,
                weight_decay=opt_cfg.weight_decay,
                betas=(opt_cfg.beta1, opt_cfg.beta2),
            )
        else:
            raise NotImplementedError(f'Optimizer not supported: {opt_cfg.type}')
        scheduler = utils_train.get_scheduler(self.config.train.scheduler, optimizer)

        if isinstance(scheduler, torch.optim.lr_scheduler.ReduceLROnPlateau):
            return {
                'optimizer': optimizer,
                'lr_scheduler': {
                    'scheduler': scheduler,
                    'monitor': 'validation/loss',
                    'frequency': int(getattr(self.config.train, 'val_freq', 1)),
                }
            }
        else:
            return [optimizer], [scheduler]

    def on_train_epoch_start(self):
        self.net_cond.eval()
        if bool(getattr(self.config.train, 'freeze_backbone', False)):
            self.model.eval()
            if self.model.cross_attn_gate is not None:
                self.model.cross_attn_gate.train()
            if not bool(getattr(self.config.train, 'freeze_residue_head', True)):
                self.model.res_inference.train()
        return super().on_train_epoch_start()

    def training_step(self, batch, batch_idx):
        data_batch = batch

        results = self.model.get_diffusion_loss(
            net_cond=self.net_cond,
            data=data_batch,
            protein_pos_apo = data_batch.protein_pos,
            protein_pos_holo=data_batch.protein_pos_holo,
            protein_v=data_batch.protein_atom_feature.float(),
            batch_protein=data_batch.protein_element_batch,
            ligand_pos=data_batch.ligand_pos,
            ligand_v=data_batch.ligand_atom_feature_full,
            batch_ligand=data_batch.ligand_element_batch,
        )
        loss, loss_ligand_pos, loss_ligand_v = results['loss'], results['loss_ligang_pos'], results['loss_v']
        loss_protein_tr, loss_protein_rot, loss_protein_chi = results['loss_protein_tr'], results['loss_protein_rot'], results['loss_protein_chi']

        self.train_iterations += 1
        # Router/gate diagnostics are cheap scalars: log every step so smoke tests
        # (which run far fewer than train_report_iter iterations) can see them.
        if self.config.model.pocket_router_mode == 'cross_attn_gate':
            for key in ('loss_gate', 'router_w_mean', 'router_w_min',
                        'router_w_max', 'router_selected_per_graph'):
                val = results.get(key, None)
                if val is not None:
                    self.log(f'train/{key}', val)
        if self.train_iterations % self.config.train.train_report_iter == 0:
            perturbed_prot_pos = results['perturbed_protein_pos']
            pred_res_tr = results['pred_res_tr']
            pred_res_rot = results['pred_res_rot']
            pred_res_chi = results['pred_res_chi']
            pred_prot_pos = apply_transforms_tensor_batch(
                protein_pos=perturbed_prot_pos,
                protein_atom_name=data_batch.protein_atom_name,
                protein_atom_to_aa_name=data_batch.protein_atom_to_aa_name,
                protein_atom_to_aa_group=data_batch.protein_atom_to_aa_group,
                protein_element_batch=data_batch.protein_element_batch,
                rotations=pred_res_rot,
                translations=pred_res_tr,
                chi_update=pred_res_chi,
                chi_mask=data_batch.protein_chi_mask,
                protein_translations_batch=data_batch.protein_translations_batch,
            )
            rmsd = torch.sqrt(torch.mean((pred_prot_pos - results['p0']) ** 2, dim=-1))
            rmsd_protein_pos = torch.mean(rmsd).item()
            self.log('iteration', self.train_iterations)
            self.log('train/loss', loss)
            self.log('train/loss_ligand_pos', loss_ligand_pos)
            self.log('train/loss_v', loss_ligand_v)
            self.log('train/loss_protein_tr', loss_protein_tr)
            self.log('train/loss_protein_rot', loss_protein_rot)
            self.log('train/loss_protein_chi', loss_protein_chi)
            self.log('train/rmsd_protein_pos', rmsd_protein_pos)
        return loss

    def on_validation_epoch_start(self):
        self.all_pred_v, self.all_true_v = [], []
        self.sum_loss, self.sum_loss_ligand_pos, self.sum_loss_v, self.sum_n = 0, 0, 0, 0
        self.sum_loss_protein_tr, self.sum_loss_protein_rot, self.sum_loss_protein_chi = 0, 0, 0
        self.sum_loss_protein_pos = 0
        self.sum_rmsd_protein_pos = 0
        return super().on_validation_epoch_start()

    def validation_step(self, batch, batch_idx):
        data_batch = batch

        batch_size = data_batch.num_graphs
        for t in np.linspace(0, self.model.num_timesteps - 1, 10).astype(int):
            time_step = torch.tensor([t] * batch_size).to(self.device)
            results = self.model.get_diffusion_loss(
                net_cond=self.net_cond,
                data=data_batch,
                protein_pos_apo=data_batch.protein_pos,
                protein_pos_holo=data_batch.protein_pos_holo,
                protein_v=data_batch.protein_atom_feature.float(),
                batch_protein=data_batch.protein_element_batch,

                ligand_pos=data_batch.ligand_pos,
                ligand_v=data_batch.ligand_atom_feature_full,
                batch_ligand=data_batch.ligand_element_batch,

                time_step=time_step
            )
            loss, loss_ligand_pos, loss_v = results['loss'], results['loss_ligang_pos'], results['loss_v']
            loss_protein_tr, loss_protein_rot, loss_protein_chi = results['loss_protein_tr'], results['loss_protein_rot'], results['loss_protein_chi']

            self.sum_loss += float(loss) * batch_size
            self.sum_loss_ligand_pos += float(loss_ligand_pos) * batch_size
            self.sum_loss_v += float(loss_v) * batch_size
            self.sum_loss_protein_tr += float(loss_protein_tr) * batch_size
            self.sum_loss_protein_rot += float(loss_protein_rot) * batch_size
            self.sum_loss_protein_chi += float(loss_protein_chi) * batch_size
            self.sum_n += batch_size
            self.all_pred_v.append(results['ligand_v_recon'].detach().cpu().numpy())
            self.all_true_v.append(data_batch.ligand_atom_feature_full.detach().cpu().numpy())

        avg_loss = self.sum_loss / self.sum_n

        return avg_loss

    def _reduce_validation_sums(self):
        values = torch.tensor(
            [
                self.sum_loss,
                self.sum_loss_ligand_pos,
                self.sum_loss_v,
                self.sum_loss_protein_tr,
                self.sum_loss_protein_rot,
                self.sum_loss_protein_chi,
                self.sum_n,
            ],
            dtype=torch.float64,
            device=self.device,
        )
        if dist.is_available() and dist.is_initialized():
            dist.all_reduce(values, op=dist.ReduceOp.SUM)
        return values

    def _gather_validation_predictions(self):
        num_classes = self.ligand_featurizer.feature_dim
        local_true = (
            np.concatenate(self.all_true_v)
            if self.all_true_v
            else np.empty((0,), dtype=np.int64)
        )
        local_pred = (
            np.concatenate(self.all_pred_v, axis=0)
            if self.all_pred_v
            else np.empty((0, num_classes), dtype=np.float32)
        )

        if not (dist.is_available() and dist.is_initialized()):
            return local_true, local_pred

        gathered = [None] * dist.get_world_size() if dist.get_rank() == 0 else None
        dist.gather_object((local_true, local_pred), gathered, dst=0)
        if dist.get_rank() != 0:
            return None, None

        global_true = np.concatenate([item[0] for item in gathered], axis=0)
        global_pred = np.concatenate([item[1] for item in gathered], axis=0)
        return global_true, global_pred

    def on_validation_epoch_end(self):
        stats = self._reduce_validation_sums()
        count = stats[6].clamp_min(1.0)
        avg_loss = stats[0] / count
        avg_loss_ligand_pos = stats[1] / count
        avg_loss_v = stats[2] / count
        avg_loss_protein_tr = stats[3] / count
        avg_loss_protein_rot = stats[4] / count
        avg_loss_protein_chi = stats[5] / count

        all_true_v, all_pred_v = self._gather_validation_predictions()
        if not (dist.is_available() and dist.is_initialized()) or dist.get_rank() == 0:
            atom_auroc = get_auroc(
                all_true_v,
                all_pred_v,
                feat_mode=self.config.data.transform.ligand_atom_mode,
            )
        else:
            atom_auroc = 0.0
        atom_auroc = torch.tensor(atom_auroc, dtype=torch.float64, device=self.device)
        if dist.is_available() and dist.is_initialized():
            dist.broadcast(atom_auroc, src=0)

        self.log('epoch', self.current_epoch)
        self.log('validation/loss', avg_loss)
        self.log('validation/loss_ligand_pos', avg_loss_ligand_pos)
        self.log('validation/loss_v', avg_loss_v * 1000)
        self.log('validation/atom_auroc', atom_auroc)
        self.log('validation/loss_protein_tr', avg_loss_protein_tr)
        self.log('validation/loss_protein_rot', avg_loss_protein_rot)
        self.log('validation/loss_protein_chi', avg_loss_protein_chi)
