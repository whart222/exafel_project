#!/bin/bash -l
#SBATCH --partition=<partition>
#SBATCH --job-name=cxid9114_index
#SBATCH --time=<walltime>
#SBATCH --nodes=<nnodes>
#SBATCH --constraint=knl
#SBATCH --image=docker:mlxd/xfel:lq79
#SBATCH --mail-type=ALL
#SBATCH --mail-user=loriordan@lbl.gov
#SBATCH --volume="/global/cscratch1/sd/mlxd/DataProc/cxid9114/cxi:/reg/d/psdm/CXI;/global/cscratch1/sd/mlxd/DataProc/cxid9114/cxi:/reg/d/psdm/cxi;/global/cscratch1/sd/mlxd/DataProc/cxid9114/regg/g:/reg/g"

# submit jobs
srun -n <nproc> shifter <srun_script>
