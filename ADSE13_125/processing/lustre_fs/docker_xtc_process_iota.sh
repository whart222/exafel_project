#!/bin/bash

source /cctbx/build/setpaths.sh

RUN=$1
TRIAL=$2
RG=$3
BASE_PATH=$4
EXP=$5
CMDMODE=$6
LIMIT=$7

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

START_XTC=$(date +"%s")

# mask and metrology files are from the current dir
IN_DIR=${BASE_PATH}/input
OUT_DIR=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out
#DATA_DIR=${BASE_PATH}/../small_xtc2/${EXP}
#DATA_DIR=/global/project/projectdirs/lcls/mona/demo18/${EXP}   # Initial path given by Mona
#DATA_DIR=/global/cscratch1/sd/psdatmgr/data/psdm/cxi/${EXP}/demo/xtc # Path given by Wilko
#DATA_DIR=/global/cscratch1/sd/asmit/iota_demo/small_xtc2/${EXP}
DATA_DIR=/global/cscratch1/sd/asmit/iota_demo/${EXP}/xtc2 # path created by me during investigation of scaling issues

export PS_CALIB_DIR=$IN_DIR
export PS_SMD_N_EVENTS=1000 #1
export PS_SMD_NODES=32 # 10

# Final set of parameters for IOTA demo 
cctbx_args="input.experiment=${EXP} input.run_num=${RUN} output.logging_dir=DevNull output.output_dir=${OUT_DIR} format.cbf.invalid_pixel_mask=/tmp/mask.pickle /tmp/params_1.phil input.xtc_dir=${DATA_DIR} input.reference_geometry=/tmp/cspad_refined_1.json"


if [ "${CMDMODE}" = "debug_timestamp" ]; then
  cctbx_args="$cctbx_args debug.event_timestamp=${ts}"
fi

# IOTA timeout of 2400 seconds used for Demo 1 on 7th Dec, 2018
cctbx_args="$cctbx_args iota.timeout_cutoff_sec=2400"

##### PLEASE NOTE #####
# ONLY PRODUCTION MODE IS SUPPORTED i.e THE LAST OPTION; PLEASE BE SURE TO USE THE RIGHT PATH TO THESE SCRIPTS
# IF USING OTHER MODES

if [ "$LIMIT" -ne 0 ]; then
    cctbx_args="$cctbx_args max_events=${LIMIT}"
fi

if [ "${CMDMODE}" = "pythonprof" ]; then
    libtbx.python -m cProfile -s tottime libtbx.python /cctbx/modules/exafel_project/ADSE13_25/xtc_process_iota_srs.py ${cctbx_args}
  
elif [ "${CMDMODE}" = "strace" ]; then
    strace -ttt -f -o $$.log libtbx.python /cctbx/modules/exafel_project/ADSE13_25/xtc_process_iota_srs.py ${cctbx_args}

elif [ "${CMDMODE}" = "debug" ]; then
    #libtbx.python ${PWD}/xtc_process_iota_srs_ps2.py ${cctbx_args}
    libtbx.python ${PWD}/xtc_process.py ${cctbx_args}
    #libtbx.python ${PWD}/mpitest.py

elif [ "${CMDMODE}" = "debug_timestamp" ]; then
    #python ${PWD}/xtc_process_iota_srs_ps2.py ${cctbx_args}
    libtbx.python ${PWD}/xtc_process.py ${cctbx_args}
  
else
    echo "Running in production mode"
    libtbx.python /tmp/xtc_process_iota_srs.py ${cctbx_args}
fi

END_XTC=$(date +"%s")
ELAPSED=$((END_XTC-START_XTC))
echo IOTA_XTC_SingleRank_TimeElapsed ${ELAPSED} ${START_XTC} ${END_XTC} ${RUN}
