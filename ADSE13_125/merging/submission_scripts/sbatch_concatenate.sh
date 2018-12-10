#!/bin/bash -l 
 
#SBATCH -C knl 
#SBATCH -N 1 
#SBATCH -q regular 
#SBATCH --exclusive 
#SBATCH -J m_concat 
#SBATCH -t 24:00:00 
#SBATCH -A m2859 

CCTBX_XFEL=/global/u2/a/asmit/cctbx.xfel
source $CCTBX_XFEL/build/setpaths.sh 

for i in `seq -f "%04g" 33 6799`; do libtbx.python $CCTBX_XFEL/modules/exafel_project/ADSE13_125/merging/concatenate_merging_dbs.py src_tag=cxic0415_${i} dest_tag=cxic0415_all src_dir=.; done
 
