#!/bin/bash

# Setup NERSC environment for building
if [ "$NERSC_HOST" == "cori" ]; then
  module swap PrgEnv-intel PrgEnv-gnu
  module swap gcc gcc/4.9.3
  #module load python/2.7-anaconda
  #source activate ${1} ##Pass conda env name as argument
fi

# Set the local bin directory on path for cmake, etc, as well as include paths
export INSTALL_DIR=$PWD/strumpack_build
export prefix=$INSTALL_DIR
export PATH=$PATH:$PWD/strumpack_build/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/strumpack_build/lib:$PWD/strumpack_build/lib64

LIB_DIR="-L$PWD/strumpack_build/lib -L$PWD/strumpack_build/lib64 -L$CONDA_PREFIX/lib "
INC_DIR="-I$PWD/strumpack_build/include -I$CONDA_PREFIX/include "

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

# STRUMPACK build
# Parameterise for Cori login/Haswell, KNL, or non-Cori linux system
pushd . > /dev/null;
cd strumpack_deps/STRUMPACK
mkdir build
cd build
if [ "$NERSC_HOST" != "cori" ]; then #Not on a Cori node; use the locally built packages for all dependencies
  echo "Standard installation"
  $INSTALL_DIR/bin/cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER=$INSTALL_DIR/bin/mpic++ -DCMAKE_C_COMPILER=$INSTALL_DIR/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=$INSTALL_DIR/bin/mpifort \
  -DSCALAPACK_LIBRARIES="$INSTALL_DIR/lib/libscalapack.a" \
  -DMETIS_INCLUDES=$INSTALL_DIR/include -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a \
  -DPARMETIS_INCLUDES=$INSTALL_DIR/include -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a \
  -DSCOTCH_INCLUDES=$INSTALL_DIR/include -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
elif [[ $(getconf _NPROCESSORS_ONLN) -ne 272 ]]; then #If not 272 cores, then we are on a login or Haswell node; use local mpi and MKL ScaLAPACK. Following the STRUMPACK Github Installer code
  echo "Default Cori login node/Haswell installation"
  cmake ../  -DBUILD_SHARED_LIBS:BOOL=ON -D_GLIBCXX_USE_CXX11_ABI=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
  -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
  -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
  -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.so;$INSTALL_DIR/lib/libscotcherr.so;$INSTALL_DIR/lib/libptscotch.so;$INSTALL_DIR/lib/libptscotcherr.so"
else #Otherwise on a KNL node; Need to adjust so that the libraries used as KNL specific
  echo "KNL installation"
  cmake ../ -D_GLIBCXX_USE_CXX11_ABI=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
   -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
   -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
   -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
   -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
   -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
fi
make -j ${proc} && make install
popd > /dev/null
