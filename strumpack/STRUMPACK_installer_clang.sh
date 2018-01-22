#!/bin/bash

# Set urls for packages
STRUMPACK=https://github.com/pghysels/STRUMPACK.git
METIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
PARMETIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
SCOTCH=https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz
TCMALLOC=https://github.com/gperftools/gperftools/releases/download/gperftools-2.6.1/gperftools-2.6.1.tar.gz

LAPACK=http://www.netlib.org/lapack/lapack-3.7.1.tgz
SCALAPACK=http://www.netlib.org/scalapack/scalapack_installer.tgz
OPENBLAS=http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
CMAKE=https://cmake.org/files/v3.9/
OPENMPI=https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz

alias mpicc='OMPI_CC=clang mpicc'
alias mpic++='OMPI_CXX=clang mpic++'

# Create directory structure
if [ ! -d ./downloads ]; then
  mkdir downloads
fi
if [ ! -d ./deps ]; then
  mkdir deps
fi
if [ ! -d ./builds ]; then
  mkdir builds
fi

pushd . >/dev/null
cd downloads

# Setup NERSC environment for building
if [ "$NERSC_HOST" == "cori" ]; then
  module swap PrgEnv-intel PrgEnv-gnu
  module swap gcc/7.1.0 gcc/6.3.0
  module load llvm/5.0.0-gnu
  #module load python/2.7-anaconda
  #source activate ${1} ##Pass conda env name as argument
fi

#################################################
# Acquire STRUMPACK-sparse and all dependencies #
#################################################
if [ ! -d ../deps/STRUMPACK ]; then
  echo "Downloading STRUMPACK"
  git clone $STRUMPACK
  cp -Rf STRUMPACK ../deps
else
  pushd . > /dev/null
  cd ../deps/STRUMPACK;
  git pull --rebase
  popd > /dev/null
fi

#Download OpenMPI
if [ ! -e OPENMPI.tar.gz ]; then
  echo "Downloading OPENMPI"
  curl -L $OPENMPI -o OPENMPI.tar.gz
  mkdir ../deps/OPENMPI
  tar xvf OPENMPI.tar.gz -C ../deps/OPENMPI --strip-components 1 
fi


# METIS v5.1.0 nested dissection
if [ ! -e METIS.tar.gz ]; then
  echo "Downloading METIS"
  curl -L $METIS -o METIS.tar.gz
  tar xvf METIS.tar.gz -C ../deps
  #Update library to 64 bit width
  if [ $(getconf LONG_BIT) == "64" ]; then
    sed -i.bak 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g' ../deps/metis-5.1.0/include/metis.h;
    rm ../deps/metis-5.1.0/include/metis.h.bak
  fi
fi

# PARMETIS v4.0.3: Parallel nested disection routines
if [ ! -e PARMETIS.tar.gz ]; then
  echo "Downloading PARMETIS"
  curl -L $PARMETIS -o PARMETIS.tar.gz
  tar xvf PARMETIS.tar.gz -C ../deps
    #Update library to 64 bit width
  if [ $(getconf LONG_BIT) == "64" ]; then
    sed -i.bak 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g' ../deps/parmetis-4.0.3/metis/include/metis.h;
    rm ../deps/parmetis-4.0.3/metis/include/metis.h.bak
  fi
fi

# SCOTCH and PT-SCOTCH v6.0.4: Matrix reordering
if [ ! -e SCOTCH.tar.gz ]; then
  echo "Downloading SCOTCH"
  curl -L $SCOTCH -o SCOTCH.tar.gz
  tar xvf SCOTCH.tar.gz -C ../deps
fi

# TCMALLOC: Faster malloc/new; Can potentially also use tbbmalloc 
# Comes with gperftools (Google performance tools)
# Not yet integrated into build
if [ ! -e TCMALLOC.tar.gz ]; then
  echo "Downloading TCMALLOC"
  curl -L $TCMALLOC -o TCMALLOC.tar.gz
  tar xvf TCMALLOC.tar.gz -C ../deps
fi

# OPENBLAS: Open source BLAS library
if [ ! -e OPENBLAS.tar.gz ]; then
  echo "Downloading OPENBLAS"
  curl -L $OPENBLAS -o OPENBLAS.tar.gz
  tar xvf OPENBLAS.tar.gz -C ../deps
fi

# CMAKE
if [ ! -e CMAKE.tar.gz ]; then
  echo "Downloading CMAKE"
  if [[ "$OSTYPE" == "linux"* ]]; then
    curl -L $CMAKE/cmake-3.9.1-Linux-x86_64.tar.gz -o CMAKE.tar.gz
    export CMAKE_DIR=$PWD/../deps/cmake-3.9.1-Linux-x86_64/
  elif [[ "$OSTYPE" == "darwin"* ]];then 
    curl -L $CMAKE/cmake-3.9.1-Darwin-x86_64.tar.gz -o CMAKE.tar.gz
    export CMAKE_DIR=$PWD/../deps/cmake-3.9.1-Darwin-x86_64/CMake.app/Contents/
  fi
  tar xvf CMAKE.tar.gz -C ../deps
  cp -Rv $CMAKE_DIR/* ../builds
fi

# SCALAPACK
if [ ! -e SCALAPACK.tar.gz ]; then
  echo "Downloading SCALAPACK"
  curl -L $SCALAPACK -o SCALAPACK.tar.gz
  tar xvf SCALAPACK.tar.gz -C ../deps
fi

popd > /dev/null
pushd . > /dev/null;

# Set the local bin directory on path for cmake, etc, as well as include paths
export INSTALL_DIR=$PWD/builds
export prefix=$INSTALL_DIR
export PATH=$PATH:$PWD/builds/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/builds/lib:$PWD/builds/lib64

LIB_DIR="-L$PWD/builds/lib -L$PWD/builds/lib64 -L$CONDA_PREFIX/lib "
INC_DIR="-I$PWD/builds/include -I$CONDA_PREFIX/include "

#Append Cray include and libs to build flags where necessary
if [ "$NERSC_HOST" == "cori" ]; then
  LIB_DIR="$LIB_DIR $(cc --cray-print-opts=libs)"
  INC_DIR="$INC_DIR $(cc --cray-print-opts=cflags)"
  
  #cc --cray-print-opts=cflags ##Get includes
  #cc --cray-print-opts=libs ##Get libs
  #cc --cray-print-opts=all ##Get both and all linker elements
fi

#Get core count/2
let proc=$(getconf _NPROCESSORS_ONLN)/2
cd deps

# Build all dependencies and dependency dependencies

#OPENMPI: Optional if to be included from modules/conda
if [ "$NERSC_HOST" != "cori" ]; then
  cd ./OPENMPI
  ./configure --prefix=$INSTALL_DIR --enable-mpi-thread-multiple
  make -j$(echo ${proc}) && make install
  cd ..
else
  echo "Nothing to do"
  #alias mpicc='cc'
  #alias mpicxx='CC'
fi

cd ./metis-5.1.0;
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc=mpicc prefix=$INSTALL_DIR
else
  MPI_CC=clang OMPI_CXX=clang++ make config cc='mpicc' prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

#Needs to have mpicc on the path; can install using conda install openmpi on Mac, or mpich/openmpi on linux
cd ./parmetis-4.0.3
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc='OMPI_CC=clang mpicc' prefix=$INSTALL_DIR
else
  OMPI_CC=clang OMPI_CXX=clang++ make config cc='mpicc' cxx='mpic++' prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

# Need to choose correct Makefile from those given; prefix needs to be passed as env variable; Path to include dir needed for MPI
# Expects mpi.h to be in /usr/include; not ideal;
cd ./scotch_6.0.4/src
cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 ./Makefile.inc
sed -i.bak 's@-O3@'"-O3 $INC_DIR"'@' ./Makefile.inc #Add CFLAGS env variable into compile path for mpi headers
sed -i.bak 's/-DSCOTCH_PTHREAD//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
sed -i.bak 's/-pthread//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
if [ "$NERSC_HOST" == "cori" ]; then
  sed -i.bak 's/gcc/OMPI_CC=clang mpicc/' ./Makefile.inc #Change default compiler
  sed -i.bak 's/mpicc/OMPI_CC=clang mpicc/' ./Makefile.inc #Change default compiler for mpi
fi
rm ./Makefile.inc.bak
OMPI_CC=clang OMPI_CXX=clang++ make scotch -j$(echo ${proc}) && OMPI_CC=clang OMPI_CXX=clang++ make ptscotch -j$(echo ${proc}) && make install
cd ../..

# Install OpenBLAS; MKL might be a good option too
if [ "$NERSC_HOST" != "cori" ]; then
  cd ./OpenBLAS-0.2.20
  OMPI_CC=clang make -j$(echo ${proc}) && make PREFIX=$INSTALL_DIR install
  cd ..
fi

if [ "$NERSC_HOST" != "cori" ]; then
  cd scalapack_installer
  python ./setup.py --downall --mpibindir=$INSTALL_DIR/bin --prefix=$INSTALL_DIR
  make -j$(echo ${proc}) && make install
  cd ..
fi

# STRUMPACK build
# Parameterise for Cori login/Haswell, KNL, or non-Cori linux system
pushd . > /dev/null;
cd STRUMPACK
mkdir build
cd build
if [ "$NERSC_HOST" != "cori" ]; then #Not on a Cori node; use the locally built packages for all dependencies
  echo "Standard installation"
  $INSTALL_DIR/bin/cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER='OMPI_CXX=clang mpic++' -DCMAKE_C_COMPILER='OMPI_CC=clang mpicc' \
  -DCMAKE_Fortran_COMPILER=$INSTALL_DIR/bin/mpifort \
  -DSCALAPACK_LIBRARIES="$INSTALL_DIR/lib/libscalapack.a" \
  -DMETIS_INCLUDES=$INSTALL_DIR/include -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a \
  -DPARMETIS_INCLUDES=$INSTALL_DIR/include -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a \
  -DSCOTCH_INCLUDES=$INSTALL_DIR/include -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
elif [[ $(getconf _NPROCESSORS_ONLN) -ne 272 ]]; then #If not 272 cores, then we are on a login or Haswell node; use local mpi and MKL ScaLAPACK. Following the STRUMPACK Github Installer code
  echo "Default Cori login node/Haswell installation"
  OMPI_CXX=clang++ OMPI_CC=clang OMPI_FC=flang cmake ../ -D_GLIBCXX_USE_CXX11_ABI=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
  -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
  -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
  -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
else #Otherwise on a KNL node; Need to adjust so that the libraries used as KNL specific
  echo "KNL installation"
  OMPI_CXX=clang++ OMPI_CC=clang OMPI_FC=flang mpifort cmake ../ -D_GLIBCXX_USE_CXX11_ABI=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
   -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpifort \
   -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
   -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
   -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
   -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
fi
OMPI_CXX=clang++ OMPI_CC=clang OMPI_FC=flang make -j ${proc} && make install
popd > /dev/null
