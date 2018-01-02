#!/bin/bash

export trial=$1
export ERR_FROM_RES=$2
export AVG_UC=$3
export INCL_NEG=$4
export UC_TOL=$5
export POSTREF=$6
export SDFAC_AUTO=$7
export SDFAC_REFINE=$8
export MPI_RANK=$9
export PROG=${10}

export TAG=${trial}
if [ -d "${TAG}" ]; then
  export TAG=${TAG}
fi
mkdir -p ${TAG}
cd ${TAG}
export TARDATA=/net/$(hostname -s)/raid1/mlxd/TAR/r0*.tar
export effective_params="d_min=2.0 \
targlob=${TARDATA} \
model=../4ngz.pdb \
backend=FS \
mysql.runtag=${trial} \
scaling.report_ML=True \
pixel_size=0.11 \
nproc=1 \
postrefinement.enable=${POSTREF} \
scaling.mtz_file=../4ngz.mtz \
scaling.mtz_column_F=f(+) \
min_corr=-1.0 \
output.prefix=${trial} \
 \
cell_rejection.unit_cell=78.95,78.95,38.12,90.0,90.0,90.0 \
raw_data.sdfac_auto=${SDFAC_AUTO} \
raw_data.sdfac_refine=${SDFAC_REFINE} \
raw_data.errors_from_sample_residuals=${ERR_FROM_RES}
raw_data.error_models.sdfac_refine.random_seed=42 \
set_average_unit_cell=${AVG_UC} \
include_negatives=${INCL_NEG} \
unit_cell_length_tolerance=${UC_TOL} \
 \
"
#if [ "$PROG" ==  "dev.cxi.mpi_merge_refltable" ]; then
#  export effective_params+="raw_data.error_models.sdfac_refine.target_function=squared"
#fi
source /net/$(hostname -s)/raid1/mlxd/InDev/miniconda/bin/activate myEnv
source /net/$(hostname -s)/raid1/mlxd/InDev/build/setpaths.sh

export PATH=/net/dials/raid1/mlxd/InDev/builds/bin:$PATH
export LD_LIBRARY_PATH=/net/dials/raid1/mlxd/InDev/builds/lib:$LD_LIBRARY_PATH

export MPICH_DBG=mpi_log
export MPICH_DBG_LEVEL=VERBOSE

mpiexec -n ${MPI_RANK} ${PROG} ${effective_params} | tee ${TAG}_merge.log
exit
#source the phenix build for single node execution:
cxi.xmerge ${effective_params}

source  ~/ml/builds/sept_sprint/myDir/phenix/phenix-installer-dev-2880-source/build/setpaths.sh
phenix.xtriage ${trial}.mtz scaling.input.xray_data.obs_labels=Iobs |tee xtriage_${trial}.out

phenix.refine /net/$(hostname -s)/raid1/mlxd/InDev/modules/exafel_project/nks/Refine_4/LM14_1colorYblyso_refine_4.pdb \
refinement.output.prefix=${trial} \
xray_data.file_name=${trial}.mtz \
xray_data.labels="Iobs" \
xray_data.r_free_flags.file_name=/net/$(hostname -s)/raid1/mlxd/InDev/modules/exafel_project/nks/merge118/reflections.mtz \
xray_data.r_free_flags.label=FreeR_flag \
main.number_of_macro_cycles=3 optimize_xyz_weight=True \
optimize_adp_weight=True nproc=20 \
ordered_solvent=True ordered_solvent.mode=every_macro_cycle \
refine.adp.individual.anisotropic="element Yb" \
refinement.input.symmetry_safety_check=warning \
refine.strategy="*individual_sites *individual_sites_real_space *rigid_body *individual_adp group_adp tls *occupancies group_anomalous" \
xray_data.force_anomalous_flag_to_be_equal_to=True \
main.wavelength=1.3853 --overwrite

libtbx.python /net/$(hostname -s)/raid1/mlxd/InDev/modules/exafel_project/nks/map_height_at_atoms.py \
${trial}_001.pdb \
${trial}.mtz \
input.xray_data.labels=Iobs \
input.symmetry_safety_check=warning \
xray_data.r_free_flags.file_name=${trial}_001.mtz \
xray_data.r_free_flags.label=R-free-flags \
selection="element Yb" |tee ${trial}_peakht.log

phenix.molprobity input.pdb.file_name=${trial}_001.pdb \
molprobity.flags.clashscore=False molprobity.flags.ramalyze=False \
molprobity.flags.omegalyze=False molprobity.flags.rotalyze=False molprobity.flags.cbetadev=False \
molprobity.flags.nqh=False molprobity.flags.rna=False molprobity.flags.restraints=False \
output.coot=False output.probe_dots=False output.prefix=${trial}_molprobity

libtbx.python /net/$(hostname -s)/raid1/mlxd/InDev/modules/exafel_project/nks/json/to_json.py ./.. ${trial}

