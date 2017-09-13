#!/bin/bash

# Set urls for packages
STRUMPACK=https://github.com/pghysels/STRUMPACK.git
METIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
PARMETIS=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
SCOTCH=https://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz
TCMALLOC=https://github.com/gperftools/gperftools/releases/download/gperftools-2.6.1/gperftools-2.6.1.tar.gz

SCALAPACK=http://www.netlib.org/scalapack/scalapack_installer.tgz
OPENBLAS=http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz
CMAKE=https://cmake.org/files/v3.9/
OPENMPI=https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz

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
if [ "$NERSC_HOST" != "cori" ]; then
  module load python/2.7-anaconda
  source activate ${1} ##Pass conda env name as argument
fi


#################################################
# Acquire STRUMPACK-sparse and all dependencies #
#################################################
if [ ! -d ../deps/STRUMPACK ]; then
  echo "Downloading STRUMPACK"
  git clone $STRUMPACK
  cp -Rf STRUMPACK ../deps
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
  cp -R $CMAKE_DIR/* ../builds
fi

# OPENMPI: Open source BLAS library
if [ ! -e OPENMPI.tar.gz ]; then
  echo "Downloading OPENMPI"
  curl -L $OPENMPI -o OPENMPI.tar.gz
  tar xvf OPENMPI.tar.gz -C ../deps
fi

# SCALAPACK: 
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

LIB_DIR="-L$PWD/builds/lib -L$PWD/builds/lib64 -L$CONDA_PREFIX/lib";
INC_DIR="-I$PWD/builds/include -I$CONDA_PREFIX/include";

#Get core count/2
let proc=$(getconf _NPROCESSORS_ONLN)/4
cd deps


# Build all dependencies and dependency dependencies

#OPENMPI: Optional if to be included from modules/conda 
if [ "$NERSC_HOST" != "cori" ]; then
  cd ./openmpi-2.1.1
  ./configure --prefix=$INSTALL_DIR --enable-mpi-thread-multiple
  make -j$(echo ${proc}) && make install
  cd ..
else
  alias mpicc='cc'
  alias mpicxx='CC'
fi

cd ./metis-5.1.0;
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc=mpicc prefix=$INSTALL_DIR
else
  make config cc=cc prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

#Needs to have mpicc on the path; can install using conda install openmpi on Mac, or mpich/openmpi on linux
cd ./parmetis-4.0.3
if [ "$NERSC_HOST" != "cori" ]; then
  make config cc=mpicc prefix=$INSTALL_DIR
else
  make config cc=cc cxx=CC prefix=$INSTALL_DIR
fi
make -j$(echo ${proc}) && make install
cd ..

# Need to choose correct Makefile from those given; prefix needs to be passed as env variable; Path to include dir needed for MPI
# Expects mpi.h to be in /usr/include; not ideal; 
cd ./scotch_6.0.4/src
cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 ./Makefile.inc

sed -i.bak 's@-O3@'"-O3 $INC_DIR"'@' ./Makefile.inc #Add CFLAGS env variable into compile path for mpi headers
#sed -i.bak 's/-DSCOTCH_PTHREAD//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
#sed -i.bak 's/-pthread//' ./Makefile.inc #Disable scotch pthreads as causes issues with MPI
rm ./Makefile.inc.bak
make scotch -j$(echo ${proc}) && make ptscotch -j$(echo ${proc}) && make install
cd ../..

# Install OpenBLAS; MKL might be a good option too
if [ "$NERSC_HOST" != "cori" ]; then
  cd ./OpenBLAS-0.2.20
  make -j$(echo ${proc}) && make PREFIX=$INSTALL_DIR install
  cd ..
fi

if [ "$NERSC_HOST" != "cori" ]; then
  cd scalapack_installer
  python ./setup.py --downall --mpibindir=$INSTALL_DIR/bin --prefix=$INSTALL_DIR
  make -j$(echo ${proc}) && make install
  cd ..
fi

if [ "$NERSC_HOST" != "cori" ]; then #Not on a Cori node; use the locally built packages for all dependencies
  echo "Standard installation"
  pushd . > /dev/null;
  cd STRUMPACK
  mkdir build 
  cd build
  $INSTALL_DIR/bin/cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER=$INSTALL_DIR/bin/mpic++ -DCMAKE_C_COMPILER=$INSTALL_DIR/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=$INSTALL_DIR/bin/mpifort \
  -DSCALAPACK_LIBRARIES="$INSTALL_DIR/lib/libscalapack.a" \
  -DMETIS_INCLUDES=$INSTALL_DIR/include -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a \
  -DPARMETIS_INCLUDES=$INSTALL_DIR/include -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a \
  -DSCOTCH_INCLUDES=$INSTALL_DIR/include -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
  make -j$(echo ${proc}) && make install
  popd > /dev/null
elif [[ $(getconf _NPROCESSORS_ONLN) -ne 272 ]]; then #If not 272 cores, then we are on a login or Haswell node; use local mpi and MKL ScaLAPACK. Following the STRUMPACK Github Installer code
  echo "KNL installation"
  pushd . > /dev/null;
  cd STRUMPACK
  mkdir build 
  cd build
  cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
  -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
  -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
  -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
  popd > /dev/null
else #Otherwise on a KNL node; Need to adjust so that the libraries used as KNL specific
  echo "Default Cori installation"
  pushd . > /dev/null;
  cd STRUMPACK
  mkdir build 
  cd build
  cmake ../ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
   -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
   -DCMAKE_EXE_LINKER_FLAGS="-dynamic" -DMETIS_INCLUDES=$INSTALL_DIR/include \
   -DMETIS_LIBRARIES=$INSTALL_DIR/lib/libmetis.a -DPARMETIS_INCLUDES=$INSTALL_DIR/include \
   -DPARMETIS_LIBRARIES=$INSTALL_DIR/lib/libparmetis.a -DSCOTCH_INCLUDES=$INSTALL_DIR/include \
   -DSCOTCH_LIBRARIES="$INSTALL_DIR/lib/libscotch.a;$INSTALL_DIR/lib/libscotcherr.a;$INSTALL_DIR/lib/libptscotch.a;$INSTALL_DIR/lib/libptscotcherr.a"
  popd > /dev/null
fi
