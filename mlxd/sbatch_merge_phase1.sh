#!/bin/bash
#SBATCH --partition=debug
#SBATCH --job-name=merge_phase1
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --constraint=haswell
#SBATCH --image=docker:mlxd/xfel:latest
#SBATCH -A lcls
#SBATCH --mail-user=loriordan@lbl.gov
#SBATCH --mail-type=ALL
#DW jobdw type=scratch capacity=60GB access_mode=striped
#DW stage_in source=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/TAR_95-114 destination=$DW_JOB_STRIPED/TAR_95-114 type=directory
#DW stage_out source=$DW_JOB_STRIPED/merge_multi destination=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/bb type=directory

#cd into the directory with the 4ngz* files
cd $1

#Set the path to the data on the burst buffer. Staging location can be changed by modifying DW stage_in above
export TARDATA=$DW_JOB_STRIPED/TAR_95-114/r0*.tar
mkdir $DW_JOB_STRIPED/merge_multi

#Copy the necessary model files to the burst buffer directory. Modify the path here to your directory with 4ngz*
cp /global/cscratch1/sd/mlxd/sept_sprint/merge_multi/4ngz* $DW_JOB_STRIPED/merge_multi
export MERGE_ROOT=$DW_JOB_STRIPED/merge_multi
export TAG=myTag
export MULTINODE=True
mkdir -p ${MERGE_ROOT}/${TAG}

srun -n 32 -c 2 --cpu_bind=cores shifter ./merge_cori.sh
