#!/bin/bash
START=$(date +"%s")

#cctbx
#source /build/setpaths.sh

#output directory
OUT_DIR=${SCRATCH}/project1/work5

libtbx.python /global/cscratch1/sd/robbol/project1/modules/cctbx_project/xfel/merging/application/input/calculate_file_load.py input.path=/global/cscratch1/sd/robbol/PSII/m*/r*/combine_experiments_*/intermediates input.reflections_suffix=reintegrated_reflections.pickle input.experiments_suffix=reintegrated_experiments.json output.output_dir=${OUT_DIR} input.parallel_file_load.method=node_memory

END=$(date +"%s")
ELAPSED=$((END-START))
echo TotalElapsed ${ELAPSED} ${START} ${END}

