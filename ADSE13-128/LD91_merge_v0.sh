#!/bin/bash

source ${CONDA_ROOT}/miniconda/bin/activate myEnv
mkdir -p ${MERGE_ROOT}/${TAG}
cd ${MERGE_ROOT}/${TAG}
export trial=${TAG}

export effective_params_nks="d_min=2.0 \
output.n_bins=10 \
model=${MERGE_ROOT}/4ngz.pdb \
targlob=${TARDATA}\
backend=FS \
pixel_size=0.11 \
nproc=64 \
postrefinement.enable=True \
postrefinement.algorithm=rs \
min_corr=-1.0 \
include_negatives=True \
set_average_unit_cell=True \
unit_cell_length_tolerance=0.02 \
raw_data.errors_from_sample_residuals=True \
lattice_rejection.unit_cell=79.0,79.0,38.3,90,90,90 \
lattice_rejection.space_group=P43212 \
scaling.mtz_file=${MERGE_ROOT}/4ngz.mtz \
scaling.mtz_column_F=f(+) \
output.prefix=${trial}"

cxi.merge ${effective_params_nks}

cxi.xmerge ${effective_params_nks}
