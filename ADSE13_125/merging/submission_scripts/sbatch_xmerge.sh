#!/bin/bash -l 
 
#SBATCH -C knl 
#SBATCH -N 1 
#SBATCH -q regular 
#SBATCH -J m_xmerge 
#SBATCH -t 24:00:00 
#SBATCH -A m2859 
 
CCTBX_XFEL=/global/u2/a/asmit/cctbx.xfel
source $CCTBX_XFEL/build/setpaths.sh
./xmerge.sh all 
