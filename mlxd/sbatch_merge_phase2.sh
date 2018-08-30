#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=merge_phase2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --constraint=haswell
#SBATCH --image=docker:mlxd/xfel:latest
#SBATCH -A lcls
#SBATCH --mail-user=loriordan@lbl.gov
#SBATCH --mail-type=ALL

cd $1

#Using SCRATCH now, as BB is not required
export MERGE_ROOT=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi
export TAG=myTag

#Copy data from the BB staged-out directory to the current working directory
cp -r ${MERGE_ROOT}/bb/${TAG}/* ${MERGE_ROOT}/${TAG}

export MULTINODE=False

#Source the Phenix ENV
source /global/project/projectdirs/lcls/mlxd/cctbx_merge/setup_env.sh
./merge_cori.sh
