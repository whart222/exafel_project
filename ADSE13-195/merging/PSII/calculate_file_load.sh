#!/bin/bash
START=$(date +"%s")

# PSII data directory
PSII_DATA_DIR=

# output directory
OUT_DIR=$PWD

export input_params="\
input.path=${PSII_DATA_DIR}/m*/r*/combine_experiments_*/intermediates \
input.reflections_suffix=reintegrated_reflections.pickle \
input.experiments_suffix=reintegrated_experiments.json \
input.parallel_file_load.method=node_memory \
input.parallel_file_load.node_memory.architecture=summit \
input.parallel_file_load.node_memory.limit=480.0 \
input.parallel_file_load.node_memory.pickle_to_memory=3.5 \
input.parallel_file_load.ranks_per_node=42 \
output.output_dir=${OUT_DIR}"

libtbx.python ${CCTBX_PREFIX}/modules/cctbx_project/xfel/merging/application/input/file_load_calculator.py ${input_params}

END=$(date +"%s")
ELAPSED=$((END-START))
echo TotalElapsed ${ELAPSED} ${START} ${END}
