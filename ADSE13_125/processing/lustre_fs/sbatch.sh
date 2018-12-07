#!/bin/bash -l

# change -N,-q, -t, NODES, RUN, TRIAL, BASE_PATH

#!/bin/bash -l
#SBATCH --account=m2859
#SBATCH --job-name=ps2cctbx
#SBATCH --nodes=100
#SBATCH --constraint=knl,quad,cache
#SBATCH --time=03:00:00
#SBATCH --image=registry.services.nersc.gov/asmit/iota_ps2_v2:latest
#SBATCH --qos=premium

NODES=100
NUM_RANKS=$((NODES*68))

t_start=`date +%s`
RUN=50 #$1
TRIAL=14 #$2
RG=2 #$3
BASE_PATH=/global/cscratch1/sd/asmit/iota_demo/cxic0415
EXP=cxic0415

echo 'Processing LCLS run ' ${RUN}

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out/
mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout/

sbcast -p ${BASE_PATH}/input/process_batch_conventional.phil /tmp/params_1.phil
sbcast -p ${BASE_PATH}/input/mask.pickle /tmp/mask.pickle
# Use refined geometry; best one from conventional indexing
sbcast -p ${BASE_PATH}/input/cspad_refined_1.json /tmp/cspad_refined_1.json
sbcast -p ${BASE_PATH}/processing/command_line/xtc_process.py /tmp/xtc_process.py

#export PMI_MMAP_SYNC_WAIT_TIME=600
t_start=`date +%s`
srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores shifter ${BASE_PATH}/processing/lustre_fs/docker_xtc_process.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} production 0
t_end=`date +%s`
echo IOTA_XTC_JobCompleted_TimeElapsed $((t_end-t_start)) $t_start $t_end ${RUN}
