#!/bin/bash -l

# change -N,-q, -t, NODES, RUN, TRIAL

#
###  !!!!  NOTE NOTE. Please read this before proceeding !!!!!! #####
# a) Note you have to make sure the .conf files write out to the right locations. It won't out of the box
# Read the NERSC documentation for more information
# b) This script stages data into and out of burst buffer.
# If you want to only stage out data, please replace the .conf file in the --bbf line with stage_out.conf
# c) You will have to have a copy of the xtc2 stream (or appropriate permissions) for the demo to be able to stage in data.

#!/bin/bash -l
#SBATCH --account=m2859
#SBATCH --job-name=iota_ps2cctbx
#SBATCH --nodes=3072
#SBATCH --constraint=knl,quad,cache
#SBATCH --time=2:00:00
#SBATCH --image=registry.services.nersc.gov/asmit/iota_ps2_v2:latest
#SBATCH --qos=regular
#SBATCH --bbf=stage_in_and_out.conf

NODES=3072
NUM_RANKS=$((NODES*68))

t_start=`date +%s`
RUN=50 #$1
TRIAL=19 #$2
RG=3 #$3
BASE_PATH=/global/cscratch1/sd/asmit/iota_demo/cxic0415
EXP=cxic0415

echo 'Processing LCLS run ' ${RUN}

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

mkdir $DW_JOB_STRIPED/out

sbcast -p ${BASE_PATH}/input/process_batch.phil /tmp/params_1.phil
sbcast -p ${BASE_PATH}/input/mask.pickle /tmp/mask.pickle
# Use refined geometry; best one from conventional indexing
sbcast -p ${BASE_PATH}/input/cspad_refined_1.json /tmp/cspad_refined_1.json
sbcast -p ${BASE_PATH}/processing/command_line/xtc_process_iota_srs.py /tmp/xtc_process_iota_srs.py

#export PMI_MMAP_SYNC_WAIT_TIME=600
t_start=`date +%s`
srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores shifter ${BASE_PATH}/processing/burst_buffer/docker_xtc_process_iota_bb.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} production 0
t_end=`date +%s`
echo IOTA_XTC_JobCompleted_TimeElapsed $((t_end-t_start)) $t_start $t_end ${RUN}
