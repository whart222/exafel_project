# Dockerfile for building cctbx.xfel with PSANA-2
# Just run 'docker build .' (please include the dot) and wait for the image to build (might take some time)
# Currently MPI libraries are being hardcoded in the LD_LIBRARY_PATH for Cori PII specific applications
FROM centos:centos7
MAINTAINER Asmit Bhowmick <abhowmick@lbl.gov>

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

    conda install -y -c lcls-ii psana && \
    conda install -y future ipython reportlab h5py mpich2 wxpython pillow libtiff mock pytest jinja2 \
      scikit-learn tabulate tqdm pytables scons mysql-python && \
    conda install -y --channel conda-forge orderedset && \
    python -m pip install procrunner && \
    rm -rf /opt/conda/lib/python2.7/site-packages/numexpr-2.6.2-py2.7.egg-info && \
    conda clean -y --all

SHELL ["/bin/bash", "-c"]
ENV BASH_ENV ~/.bashrc
ENV PYTHONPATH=/opt/conda/lib/python2.7/site-packages
# build cctbx
RUN mkdir cctbx &&\
        cd cctbx &&\
        wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py &&\
        python bootstrap.py hot update build --builder=dials --with-python=/opt/conda/bin/python --nproc=4 &&\
        cd .. &&\
        conda uninstall --force mpich2 && \
        yum -y install strace

# setup environments for psana2 and cctbx
ENV CONDA_EXE=/opt/conda/bin/conda
ENV CONDA_PREFIX=/opt/conda
ENV CONDA_PROMPT_MODIFIER=(base)
ENV CONDA_PYTHON_EXE=/opt/conda/bin/python
ENV CONDA_SHLVL=1
ENV INFOPATH=/opt/rh/devtoolset-7/root/usr/share/info

ENV LD_LIBRARY_PATH=/opt/udiImage/modules/mpich/lib64:/opt/rh/devtoolset-7/root/usr/lib64:/opt/rh/devtoolset-7/root/usr/lib:/opt/rh/devtoolset-7/root/usr/lib64/dyninst:/opt/rh/devtoolset-7/root/usr/lib/dyninst:/opt/rh/devtoolset-7/root/usr/lib64:/opt/rh/devtoolset-7/root/usr/lib
ENV MANPATH=/opt/rh/devtoolset-7/root/usr/share/man:/usr/common/software/man:/usr/common/mss/man:/usr/common/nsg/man:/opt/cray/pe/mpt/7.7.0/gni/man/mpich:/opt/cray/pe/atp/2.1.1/man:/opt/cray/alps/6.5.28-6.0.5.0_18.6__g13a91b6.ari/man:/opt/cray/job/2.2.2-6.0.5.0_8.47__g3c644b5.ari/man:/opt/cray/pe/pmi/5.0.13/man:/opt/cray/pe/libsci/18.03.1/man:/opt/cray/pe/man/csmlversion:/opt/cray/pe/craype/2.5.14/man:/opt/intel/compilers_and_libraries_2018.1.163/linux/man/common:/usr/syscom/nsg/man:/opt/cray/pe/modules/3.2.10.6/share/man:/global/homes/c/canon/man:/usr/local/man:/usr/share/man:/opt/cray/share/man:/opt/cray/pe/man:/opt/cray/share/man
ENV PATH=/cctbx/build/bin:/lcls2/install/bin:/lcls2/build/bin:/opt/rh/devtoolset-7/root/usr/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/udiImage/bin
ENV PCP_DIR=/opt/rh/devtoolset-7/root
ENV PERL5LIB=/opt/rh/devtoolset-7/root//usr/lib64/perl5/vendor_perl:/opt/rh/devtoolset-7/root/usr/lib/perl5:/opt/rh/devtoolset-7/root//usr/share/perl5/vendor_perl
ENV PMI_MMAP_SYNC_WAIT_TIME=600
ENV PYTHONPATH=/opt/conda/lib/python2.7/site-packages
ENV HOME=/tmp

# Clone exafel_project and configure it
RUN cd /cctbx/modules/ && \
  git clone https://github.com/ExaFEL/exafel_project.git &&\
  cd /cctbx/build/ && source setpaths.sh && libtbx.configure exafel_project && make
