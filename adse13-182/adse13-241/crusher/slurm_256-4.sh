#!/bin/bash -l

#SBATCH -p batch            # partition
#SBATCH -N 256                # Number of nodes
#SBATCH -t 00:15:00         # wall clock time limit
#SBATCH -J fom_kokkos_256      # job name
#SBATCH -A chm137_crusher   # allocation
#SBATCH --gpus-per-node 4   # devices per node
#SBATCH -o job%j.out
#SBATCH -e job%j.err
# do not use #SBATCH --exclusive

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export N_SIM=100000 # total number of images to simulate
export ADD_BACKGROUND_ALGORITHM=cuda # cuda or jh or sort_stable
export DEVICES_PER_NODE=4
export MOS_DOM=25

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd;ls
srun -n 4096 -c 16 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/revapi/step5_batch.py context=kokkos_gpu
echo "jobend $(date)";pwd;ls
