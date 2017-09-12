#!/bin/bash -l
#SBATCH --partition=regular
#SBATCH --job-name=cxid9114_mask
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --image=docker:mlxd/xfel:lq79
#SBATCH --mail-type=ALL
#SBATCH --mail-user=loriordan@lbl.gov
#SBATCH --volume="/global/cscratch1/sd/mlxd/DataProc/cxid9114/cxi:/reg/d/psdm/CXI;/global/cscratch1/sd/mlxd/DataProc/cxid9114/cxi:/reg/d/psdm/cxi;/global/cscratch1/sd/mlxd/DataProc/cxid9114/regg/g:/reg/g"

if [ ! -f ./mask.sh ]; then
    echo "#\!/bin/bash" >> ./mask.sh
    echo "export SIT_DATA=/reg/g/psdm/data"  >> ./mask.sh
    echo "export SIT_PSDM_DATA=/reg/d/psdm"  >> ./mask.sh
    echo "source /build/setpaths.sh"  >> ./mask.sh
    echo "cxi.make_dials_mask --maxproj-min=50 -o mask.pickle \${1} \${2} \${3}" >> ./mask.sh
fi

# $1 = Dark average, $2 = Dark stddev, $3 = bright max
srun -n 68 -c 1 --cpu_bind=cores shifter ./mask.sh ${1} ${2} ${3}
