#!/bin/bash

# Set urls for packages
STRUMPACK=https://github.com/pghysels/STRUMPACK.git
METIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
PARMETIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
SCOTCH=https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz
TCMALLOC=https://github.com/gperftools/gperftools/releases/download/gperftools-2.6.1/gperftools-2.6.1.tar.gz

LAPACK=http://www.netlib.org/lapack/lapack-3.8.0.tar.gz
SCALAPACK=http://www.netlib.org/scalapack/scalapack_installer.tgz
OPENBLAS=http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
CMAKE=https://cmake.org/files/v3.9/
OPENMPI=https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz

# Create directory structure
if [ ! -d ./strumpack_downloads ]; then
  mkdir strumpack_downloads
fi
if [ ! -d ./strumpack_deps ]; then
  mkdir strumpack_deps
fi
if [ ! -d ./strumpack_build ]; then
  mkdir strumpack_build
fi

pushd . >/dev/null
cd strumpack_downloads

# Setup NERSC environment for building
if [ "$NERSC_HOST" == "cori" ]; then
  module swap PrgEnv-intel PrgEnv-gnu
  module swap gcc gcc/4.9.3
  #module load python/2.7-anaconda
  #source activate ${1} ##Pass conda env name as argument
fi

#################################################
# Acquire STRUMPACK-sparse and all dependencies #
#################################################
if [ ! -d ../strumpack_deps/STRUMPACK ]; then
  echo "Downloading STRUMPACK"
  git clone $STRUMPACK
  cp -Rf STRUMPACK ../strumpack_deps
else
  pushd . > /dev/null
  cd ../strumpack_deps/STRUMPACK;
  git pull --rebase
  popd > /dev/null
fi

#Download OpenMPI
if [ ! -e OPENMPI.tar.gz ]; then
  echo "Downloading OPENMPI"
  curl -L $OPENMPI -o OPENMPI.tar.gz
  mkdir ../strumpack_deps/OPENMPI
  tar xvf OPENMPI.tar.gz -C ../strumpack_deps/OPENMPI --strip-components 1 
fi


# METIS v5.1.0 nested dissection
if [ ! -e METIS.tar.gz ]; then
  echo "Downloading METIS"
  curl -L $METIS -o METIS.tar.gz
  tar xvf METIS.tar.gz -C ../strumpack_deps
  #Update library to 64 bit width
  if [ $(getconf LONG_BIT) == "64" ]; then
    sed -i.bak 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g' ../strumpack_deps/metis-5.1.0/include/metis.h;
    rm ../strumpack_deps/metis-5.1.0/include/metis.h.bak
  fi
fi

# PARMETIS v4.0.3: Parallel nested disection routines
if [ ! -e PARMETIS.tar.gz ]; then
  echo "Downloading PARMETIS"
  curl -L $PARMETIS -o PARMETIS.tar.gz
  tar xvf PARMETIS.tar.gz -C ../strumpack_deps
    #Update library to 64 bit width
  if [ $(getconf LONG_BIT) == "64" ]; then
    sed -i.bak 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g' ../strumpack_deps/parmetis-4.0.3/metis/include/metis.h;
    rm ../strumpack_deps/parmetis-4.0.3/metis/include/metis.h.bak
  fi
fi

# SCOTCH and PT-SCOTCH v6.0.4: Matrix reordering
if [ ! -e SCOTCH.tar.gz ]; then
  echo "Downloading SCOTCH"
  curl -L $SCOTCH -o SCOTCH.tar.gz
  tar xvf SCOTCH.tar.gz -C ../strumpack_deps
fi

# TCMALLOC: Faster malloc/new; Can potentially also use tbbmalloc 
# Comes with gperftools (Google performance tools)
# Not yet integrated into build
if [ ! -e TCMALLOC.tar.gz ]; then
  echo "Downloading TCMALLOC"
  curl -L $TCMALLOC -o TCMALLOC.tar.gz
  tar xvf TCMALLOC.tar.gz -C ../strumpack_deps
fi

# OPENBLAS: Open source BLAS library
if [ ! -e OPENBLAS.tar.gz ]; then
  echo "Downloading OPENBLAS"
  curl -L $OPENBLAS -o OPENBLAS.tar.gz
  tar xvf OPENBLAS.tar.gz -C ../strumpack_deps
fi

# CMAKE
if [ ! -e CMAKE.tar.gz ]; then
  echo "Downloading CMAKE"
  if [[ "$OSTYPE" == "linux"* ]]; then
    curl -L $CMAKE/cmake-3.9.1-Linux-x86_64.tar.gz -o CMAKE.tar.gz
    export CMAKE_DIR=$PWD/../strumpack_deps/cmake-3.9.1-Linux-x86_64/
  elif [[ "$OSTYPE" == "darwin"* ]];then 
    curl -L $CMAKE/cmake-3.9.1-Darwin-x86_64.tar.gz -o CMAKE.tar.gz
    export CMAKE_DIR=$PWD/../strumpack_deps/cmake-3.9.1-Darwin-x86_64/CMake.app/Contents/
  fi
  tar xvf CMAKE.tar.gz -C ../strumpack_deps
  cp -Rv $CMAKE_DIR/* ../strumpack_build
fi

# SCALAPACK
if [ ! -e SCALAPACK.tar.gz ]; then
  echo "Downloading SCALAPACK"
  curl -L $SCALAPACK -o SCALAPACK.tar.gz
  tar xvf SCALAPACK.tar.gz -C ../strumpack_deps
fi

popd > /dev/null
pushd . > /dev/null;

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
cd strumpack_deps

# Build all dependencies and dependency dependencies

#OPENMPI: Optional if to be included from modules/conda
cd ./OPENMPI
./configure --prefix=$INSTALL_DIR --enable-static --disable-shared --with-verbs=no CFLAGS='-fPIC -m64' CXXFLAGS='-fPIC -m64' FCFLAGS='-fPIC -m64' --disable-oshmem #Needed to build on Cray systems
make -j$(echo ${proc}) && make install
cd ..

cd ./metis-5.1.0;
#Modify the shared flag to emit position independent code only. Shared library not required, but fPIC static lib is necessary.
sed -i.bak 's/-DSHARED=1/-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true/' ./Makefile
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc=mpicc shared=1 prefix=$INSTALL_DIR
else
  make config cc=cc shared=1 prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

#Needs to have mpicc on the path; can install using conda install openmpi on Mac, or mpich/openmpi on linux
cd ./parmetis-4.0.3
sed -i.bak 's/-DSHARED=1/-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true/' ./Makefile
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc=mpicc shared=1 prefix=$INSTALL_DIR
else
  make config cc=cc shared=1 cxx=CC prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

# Need to choose correct Makefile from those given; prefix needs to be passed as env variable; Path to include dir needed for MPI
# Expects mpi.h to be in /usr/include; not ideal;
cd ./scotch_6.0.4/src
cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 ./Makefile.inc
sed -i.bak 's@-O3@'"-O3 -fPIC $INC_DIR"'@' ./Makefile.inc #Add CFLAGS env variable into compile path for mpi headers
sed -i.bak 's/-DSCOTCH_PTHREAD//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
sed -i.bak 's/-pthread//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
if [ "$NERSC_HOST" == "cori" ]; then
  sed -i.bak 's/gcc/cc/' ./Makefile.inc #Change default compiler
  sed -i.bak 's/mpicc/cc/' ./Makefile.inc #Change default compiler for mpi
fi
rm ./Makefile.inc.bak
make scotch -j$(echo ${proc}) && make ptscotch -j$(echo ${proc}) && make prefix=${INSTALL_DIR} install
cd ../..

# Install OpenBLAS; MKL might be a good option too
#if [ "$NERSC_HOST" != "cori" ]; then
  cd ./OpenBLAS-0.2.20
  make -j$(echo ${proc}) FC='ftn -fPIC' CC='cc -fPIC' NO_SHARED=1 && make PREFIX=$INSTALL_DIR install
  cd ..
#fi

#if [ "$NERSC_HOST" != "cori" ]; then
  cd scalapack_installer
  python ./setup.py --downall --ccflags='-O3 -fPIC' --mpibindir=$INSTALL_DIR/bin --prefix=$INSTALL_DIR
  make -j$(echo ${proc}) && make install
  cd ..
#fi

# STRUMPACK build
# Parameterise for Cori login/Haswell, KNL, or non-Cori linux system
pushd . > /dev/null;
cd STRUMPACK
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
  cmake ../ -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -D_GLIBCXX_USE_CXX11_ABI=0 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
  -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
  -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
  -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
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
