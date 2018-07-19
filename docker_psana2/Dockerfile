FROM centos:centos7

RUN yum clean all && \
    yum -y update && \
    yum -y install centos-release-scl-rh && \
    yum -y install devtoolset-7-gcc devtoolset-7-gcc-c++ devtoolset-7-gcc-gfortran \
                   devtoolset-7-make wget bzip2 && \
    yum clean all -y

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    echo "source /opt/rh/devtoolset-7/enable" >> ~/.bashrc && \
    source ~/.bashrc && \
    conda update -y conda && \
    conda update -y --all && \
    conda install numpy cython cmake git && \
    conda install mpi4py future h5py mpich2 wxpython pillow \
                  libtiff mock pytest jinja2 scikit-learn tabulate &&\
    conda install --channel conda-forge orderedset &&\
    python -m pip install procrunner &&\
    conda clean -y --all

SHELL ["/bin/bash", "-c"]
ENV BASH_ENV ~/.bashrc

# build psana2 with python 2
RUN git clone https://github.com/slac-lcls/lcls2.git && \
    cd lcls2 && \
    ./build_python2_psana.sh

RUN mkdir -p /reg/g &&\
    mkdir -p /reg/d/psdm/CXI &&\
    mkdir -p /reg/d/psdm/cxi &&\
    mkdir -p /reg/d/psdm/mfx &&\
    mkdir -p /reg/d/psdm/MFX

# cctbx
RUN mkdir cctbx &&\
    cd cctbx &&\
    wget https://raw.githubusercontent.com/monarin/cctbx_project/master/libtbx/auto_build/bootstrap.py &&\
    python bootstrap.py hot update build --builder=dials --with-python=/opt/conda/bin/python --nproc=1 &&\
    conda uninstall --force mpich2 &&\
    cd ..
