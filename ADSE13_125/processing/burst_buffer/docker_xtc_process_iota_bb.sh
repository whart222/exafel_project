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
OUT_DIR=${DW_JOB_STRIPED}/out
DATA_DIR=$DW_JOB_STRIPED # When xtc2 streams are staged to Burst Buffer

export PS_CALIB_DIR=$IN_DIR
export PS_SMD_N_EVENTS=1000 #1
export PS_SMD_NODES=32 # 10

# xtc_process basic params to be specified here
cctbx_args="input.experiment=${EXP} input.run_num=${RUN} output.logging_dir=DevNull output.output_dir=${OUT_DIR} format.cbf.invalid_pixel_mask=/tmp/mask.pickle /tmp/params_1.phil input.xtc_dir=${DATA_DIR} input.reference_geometry=/tmp/cspad_refined_1.json"


if [ "${CMDMODE}" = "debug_timestamp" ]; then
  cctbx_args="$cctbx_args debug.event_timestamp=${ts}"
fi

# Timeout of 2400 seconds used for Demo 1 
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
