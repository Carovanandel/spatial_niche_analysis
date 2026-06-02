
# Sparrow env installation
To create the sparrow environment, use the included sparrow_env.yml file, which includes the link to the exact sparrow git commit needed. \
Installation with pip install sparrow does not work, as sparrow is not on PyPi. \
Numpy v.1.26.4 is needed, but this gives some dependency issues. Installing it in this way (no specific version included in the .yml, then downgrading to 1.26.4 after installing the .yml) gives a pip dependency warning, but creates a working environment. \

```
conda create -n sparrow_copy python=3.10 pip
chmod -R u+w conda/envs/sparrow_copy
conda activate sparrow_copy
conda env update -n sparrow_copy -f sparrow_env.yml --prune
pip install --upgrade numpy==1.26.4
conda install ipykernel
```

# CellCharter env installation
```
module load tools/miniconda/python3.10/23.3.1 
conda activate
conda create -n cellcharter -c conda-forge python=3.10 mamba
chmod -R u+w conda/envs/cellcharter
conda activate cellcharter
mamba install pytorch torchvision torchaudio pytorch-cuda=12.4 -c pytorch -c nvidia -c conda-forge
pip install scvi-tools
pip install cellcharter
conda install ipykernel
conda install -c conda-forge python-igraph
conda install -c conda-forge leidenalg
pip install pygam      ### OPTIONAL: package for using generalized additive models
```

# scArches/trVAE env installation
```
conda create -n scarches python=3.10 pip
chmod -R u+w /exports/archive/hg-groep-spitali/Students/Caro/conda/envs/scarches
conda activate scarches
conda install pytorch=2.0.1 torchvision=0.15.2 torchaudio=2.0.2 pytorch-cuda=11.7 numpy=1.26.4 -c pytorch -c nvidia
pip install scvi-tools scarches
conda install ipykernel
```

# Palantir env installation
```
conda create -n palantir python=3.10 pip
chmod -R u+w /exports/archive/hg-groep-spitali/Students/Caro/conda/envs/palantir
conda activate palantir
pip install palantir
conda install ipykernel
```