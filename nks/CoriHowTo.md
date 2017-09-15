
# Cori install/usage instructions

Installation of the Conda environment 
-------------------------------------
Instructions for 0917-sprint-CCTBX + SCIPY build are based on https://exafel.github.io/docs/psana-cctbx-install.

1. Create and change to a new directory in bash for installation of Conda and CCTBX, setting the environment variable $CONDA_ROOT to this path.
```bash
mkdir myDir; 
export CONDA_ROOT=$PWD/myDir; 
cd $CONDA_ROOT
```
2. Install Miniconda for Python 2.7 (64-bit) to directory `$CONDA_ROOT/miniconda`, acquire packages and create environment for building cctbx/Phenix [do not prepend to .bashrc]
```bash
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh; 
./Miniconda2-latest-Linux-x86_64.sh -b -p $CONDA_ROOT/miniconda;
$CONDA_ROOT/miniconda/bin/conda update -y conda;
$CONDA_ROOT/miniconda/bin/conda create -n myEnv -y;
source $CONDA_ROOT/miniconda/bin/activate myEnv;
$CONDA_ROOT/miniconda/bin/conda install -y --channel lcls-rhel7 psana-conda;
$CONDA_ROOT/miniconda/bin/conda install -y IPython h5py mpich2 wxpython pillow libtiff mysql-python
```

Installation of CCTBX-XFEL
--------------------------
3. Change to the $CONDA_ROOT directory, verify the Conda Python is default, and acquire the bootstrap file for building CCTBX in the sept_sprint branch:
```bash 
cd $CONDA_ROOT
which python
curl -L https://raw.githubusercontent.com/ExaFEL/cctbx_project/sept_sprint/libtbx/auto_build/bootstrap.py -o bootstrap.py
```

4. Pull all the required dependencies for building CCTBX into the `modules` directory, then build into `build`.
```bash
python bootstrap.py hot update --builder=xfel --sfuser=<github user name>
python bootstrap.py build --builder=xfel --with-python=`which python` --nproc=<# cores available; 8 for Cori login nodes is fine> 
```
5. Source the newly built CCTBX paths [Note: Conda MUST be sourced path before running the following command. This applies when opening a new tmerinal session]
```bash
source $CONDA_ROOT/build/setpaths.sh;
export BOOST_ADAPTBX_FPE_DEFAULT=1 # When using matplotlib, always set this. Ignore otherwise.
```
Installation of separate PHENIX build
-------------------------------------
6. Create a new directory for the Phenix installation (can be within $CONDA_ROOT)
```bash
cd $CONDA_ROOT
mkdir phenix; export PHENIX_ROOT=$PWD/phenix
cd $PHENIX_ROOT
```
7.  Request a download password for Phenix from http://phenix-online.org and download the tar file with these credentials:
```bash
#Ask password flag prompts for the password before downloading
wget --user=<username> --ask-password \
"https://www.phenix-online.org/download/phenix/nightly/send_octet_stream.cgi?version=dev-2880&file=phenix-installer-dev-2880-source.tar.gz";

#Filename can be mangled, so a rename is provided for simplicity
mv ./send_octet_stream.cgi\?version\=dev-2880\&file\=phenix-installer-dev-2880-source.tar.gz phenix.tar.gz
```

8. Untar the file, cd into the directory and build the binaries:
```bash
cd phenix-installer-dev-2880-source
cp modules/cctbx_project/libtbx/auto_build/bootstrap.py .;
cat ./bootstrap.py | sed -e '1515,1529d;' > bootstrap_2.py ; #Remove gui error
mv bootstrap_2.py ./bootstrap.py;
python bootstrap.py build --builder=phenix --with-python=`which python` --nproc=<cores available for compile> 
```

Running the LD91, run 108 example
---------------------------------
A. We must now set environment variables for the database access, and creating the appropriate output from the merging steps:
```bash
# a unique tag for this merging trial with a particular $TARDATA.
export TAG=<any tag> 

# glob describing path to TAR files (run 108 in this case). 
# This path is currently pointing to the GPFS filesystem, 
# so copying to SCRATCH is recommended for optimal performance.
export TARDATA=/global/project/projectdirs/lcls/mlxd/cxid9114/processing/batch_metrology_r0113_013/TAR_95-114/r0108*.tar 

# Path for the top-level merging directory above several $TAG/$TARDATA pairs. 
# Preferably somewhere on SCRATCH
export MERGE_ROOT=<new directory name> 

# This file will set all the required database parameters.
# A path to this will be provided if requested, as it will not be checked into the repository.
source db_params.dat
```

B. Next, we must next download the reference data for lysozyme from the PDB
```bash
cd $MERGE_ROOT; 
phenix.fetch_pdb --mtz 4ngz; 
cd -
```
C. The appropriate single or multi-node sections of the code are now called. We perform this interactively here for simplicity. Batch submission will be provided later. First, the multinode section is run with
```bash
# True=multinode (aka merge), False=single node (xmerge,etc). 
#Note: the script needs to be run twice. Once with MULTINODE=True for postrefinement.
#Second time with MULTINODE=False for merging stats, xtriage, phenix.refine, and anomalous stats
export MULTINODE=True #This value chooses which mode to operate within. `

#Interactively acquire a KNL node and perform the merging step with MULTINODE=True
salloc -N 1 -C knl -A m2859 -p debug -t 00:30:00 --qos=interactive --image=docker:mlxd/xfel:latest;
$CONDA_ROOT/modules/exafel_project/nks/merge_v02.sh 260 #The value of `nproc` is passed in as an argument.
exit; #Give up the interactive KNL session
```

Next, we use the single node section on Haswell. 
```bash
salloc -N 1 -C haswell -A m2859 -p debug -t 00:30:00 --qos=interactive 
export MULTINODE=False
$CONDA_ROOT/modules/exafel_project/nks/merge_v02.sh 32
exit #Give up the interactive Haswell session
```

The resulting merged data should now be placed into the $MERGE_ROOT directory.
