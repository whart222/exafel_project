#!/usr/bin/bash

# setup compiler
module load gcc/7.4.0

# install conda if it is not available
CONDA_LOCATION=${HOME}/software/miniconda3
if [ ! -d ${CONDA_LOCATION} ]; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-ppc64le.sh
  bash Miniconda3-latest-Linux-ppc64le.sh -b -p ${CONDA_LOCATION}
  rm Miniconda3-latest-Linux-ppc64le.sh
fi
source ${CONDA_LOCATION}/etc/profile.d/conda.sh

# install DIALS
DIALS_LOCATION=${HOME}/software/dials

if [ ! -d ${DIALS_LOCATION} ]; then
  mkdir -p ${DIALS_LOCATION}
fi

cd ${DIALS_LOCATION}

# get bootstrap.py
wget "https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py"

# download sources
python bootstrap.py hot update --builder=dials

# create conda_base
conda create -y -p ./conda_base python=2 biopython cython future hdf5 h5py ipython jinja2 matplotlib mock msgpack-python numpy pillow pip psutil pyopengl pytest pytest-mock pytest-xdist pyyaml reportlab scikit-learn six tabulate tqdm=4.23.4
conda activate ./conda_base
pip install dials-data
pip install mrcfile
pip install orderedset
pip install procrunner
conda deactivate

# build
python bootstrap.py build --builder=dials --use-conda --nproc=16

cd -

