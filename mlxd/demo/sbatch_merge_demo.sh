#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=merge_phase1
#SBATCH --time=00:25:00
#SBATCH --nodes=2
#SBATCH --constraint=haswell
#SBATCH --image=docker:mlxd/xfel:alpha
#SBATCH -A lcls
#SBATCH -o 2node_cs_all_null.log
#SBATCH -e 2node_cs_all_null.err
#SBATCH --mail-user=loriordan@lbl.gov
#SBATCH --mail-type=ALL
#DW jobdw type=scratch capacity=60GB access_mode=striped
#DW stage_in source=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/TAR_95-114 destination=$DW_JOB_STRIPED/TAR_95-114 type=directory
#DW stage_out source=$DW_JOB_STRIPED/merge_multi destination=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/bb type=directory

#cd into the directory with the 4ngz* files
cd $1

export TARDATA=$DW_JOB_STRIPED/TAR_95-114/r0*.tar
mkdir $DW_JOB_STRIPED/merge_multi

#Copy the necessary model files to the burst buffer directory. Modify the path here to your directory with 4ngz*
cp /global/cscratch1/sd/mlxd/sept_sprint/merge_multi/4ngz* $DW_JOB_STRIPED/merge_multi
export MERGE_ROOT=$DW_JOB_STRIPED/merge_multi
export TAG=2node_cs_all_null
export MULTINODE=True
mkdir -p ${MERGE_ROOT}/${TAG}
echo "SRUN START RANK=0 TIME=" $(date +%s)
srun -n 64 -c 2 --cpu_bind=cores shifter ${1}/merge_cori_demo.sh
echo "SRUN STOP RANK=0 TIME=" $(date +%s)

#Stage 2
export MULTINODE=False

#Source the Phenix ENV
source /global/project/projectdirs/lcls/mlxd/cctbx_merge/setup_env.sh
./merge_cori_demo.sh
