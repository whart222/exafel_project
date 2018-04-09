#!/bin/bash
. <SPACK_ROOT>/share/spack/setup-env.sh

spack install miniconda2

#Load the spack environment packages
for mods in miniconda2 mpi openblas scotch parmetis metis netlib-scalapack strumpack;
do
  echo "$(spack module loads ${mods} | grep -v '#')"
done

#Setup the conda environment
conda update -y conda;
conda create -n myEnv -y conda;
source activate myEnv
conda install -y IPython h5py wxpython pillow libtiff mysql-python jinja2 matplotlib

#Build packages
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
python bootstrap.py hot update --builder=xfel --sfuser=<USERNAME> --cciuser=<USERNAME>
python bootstrap.py build --builder=xfel --with-python=`which python` --nproc=<NUM_CORES>
source ./build/setpaths.sh

# Install mpi4py with the available MPI environment
wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
tar xvf mpi4py-3.0.0.tar.gz
cd mpi4py-3.0.0
libtbx.python setup.py build --mpicc=$(which mpicc)
libtbx.python setup.py build_exe --mpicc=$(which cc)
libtbx.python setup.py install
libtbx.python setup.py install_exe

#Create strumpack_build folder and populate with the contents of the spack bin/include/lib folders
cd ..
mkdir -p strumpack_build && cd strumpack_build
for mods in mpi openblas scotch parmetis metis netlib-scalapack strumpack;
do
  ln -sf $(spack location --install-dir ${mods})/* .
done

