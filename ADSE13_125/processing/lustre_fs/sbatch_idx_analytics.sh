#!/bin/bash -l

#SBATCH -C knl
#SBATCH -N 8
#SBATCH -q premium
#SBATCH -J idx_016
#SBATCH -t 1:00:00 
#SBATCH -A m2859 

NODES=8
NUM_RANKS=$((NODES*68))

CCTBX_XFEL=/global/u2/a/asmit/cctbx.xfel
source $CCTBX_XFEL/build/setpaths.sh

srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores libtbx.python $CCTBX_XFEL/modules/exafel_project/ADSE13_25/command_line/indexing_analytics.py mpi=True write_out_timings=False wall_time=4534 num_nodes=3500

