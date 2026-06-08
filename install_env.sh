#!/bin/bash

set -euo pipefail

ENV_NAME=${ENV_NAME:-apo2mol}

source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    conda activate "${ENV_NAME}"
else
    conda create -y -n "${ENV_NAME}" python=3.10
    conda activate "${ENV_NAME}"
fi

python -m pip install --upgrade pip setuptools wheel

pip install pytorch-lightning==2.2.3
pip install torch==2.4.0 torchvision==0.19.0 torchaudio==2.4.0 --index-url https://download.pytorch.org/whl/cu124
pip install numpy==1.24.1
pip install torch-geometric==2.0.4
pip install ninja cmake
PYG_WHEEL_URL=${PYG_WHEEL_URL:-https://data.pyg.org/whl/torch-2.4.0+cu124.html}
pip install --only-binary=:all: torch-scatter==2.1.2 torch-sparse==0.6.18 torch-cluster==1.6.3 torch-spline-conv==1.2.2 -f "${PYG_WHEEL_URL}"
pip install bypy==1.8.5
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
pip install einops==0.8.0 easydict==1.13
pip install huggingface-hub==0.25.2 hydra-core==1.3.2 ipython==8.28.0 jinja2==3.1.3
conda install -y -c conda-forge patch
pip install lmdb==1.5.1
pip install matplotlib==3.9.2 meeko==0.1.dev3 mmcif-pdbx==2.0.1
pip install omegaconf==2.3.0 pandas==2.2.3 pdb2pqr==3.6.1 pillow==10.2.0 py3dmol==2.4.0
pip install tornado==6.4.1 tqdm==4.66.5 wandb==0.18.3
conda install -y -c conda-forge vina=1.2.5
pip install wheel==0.37.1
pip install easydict==1.13
pip install biopython
pip install kornia
pip install numpy-quaternion==2022.4.4
conda install -y -c conda-forge rdkit=2022.09.5 openbabel=3.1.1
