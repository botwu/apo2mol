import torch
from torch.utils.data import Subset
import pickle
from tqdm import tqdm
from .mySubset import GetSubset
from .pl_pair_dataset import PocketLigandPairDataset


def _filter_missing_lmdb_indices(dataset, indices, split_name):
    filtered = []
    missing = []
    for idx in indices:
        idx = int(idx)
        if dataset.has_key(idx):
            filtered.append(idx)
        else:
            missing.append(idx)
    if missing:
        print(
            f'[Apo2Mol] filtered {len(missing)} missing LMDB entr'
            f'{"y" if len(missing) == 1 else "ies"} from {split_name} split; '
            f'examples: {missing[:10]}'
        )
    return filtered


def get_dataset(config, *args, **kwargs):
    name = config.name
    root = config.path
    index_path = config.index
    pocket_type = config.type

    assert 'split' in config
    split_indices_dict = torch.load(config.split, weights_only=False)
    train_split_indices = split_indices_dict['train']
    valid_split_indices = split_indices_dict['valid']
    test_split_indices = split_indices_dict['test']  # the 'key' of val_dataset is 'test'

    assert name == 'pl'
    dataset = PocketLigandPairDataset(root, index_path, pocket_type, *args, **kwargs)
    train_split_indices = _filter_missing_lmdb_indices(dataset, train_split_indices, 'train')
    valid_split_indices = _filter_missing_lmdb_indices(dataset, valid_split_indices, 'valid')
    test_split_indices = _filter_missing_lmdb_indices(dataset, test_split_indices, 'test')

    train_dataset = GetSubset(dataset, indices=train_split_indices)
    valid_dataset = GetSubset(dataset, indices=valid_split_indices)
    test_dataset = GetSubset(dataset, indices=test_split_indices)

    subsets = {'train': train_dataset, 'valid': valid_dataset, 'test': test_dataset}
    return subsets
