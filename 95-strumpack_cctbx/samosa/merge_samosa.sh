#!/bin/bash
export TAG=step4K_samosa_debug_10k_test_cc
export trial=${TAG}
export MERGE_ROOT=/net/dials/raid1/mlxd/InDev/STRUMPACK_MPI_DIST/sam/merge #/net/dials/raid1/sauter/LS49_merge
export SAMOSA_ROOT=/net/dials/raid1/mlxd/InDev/STRUMPACK_MPI_DIST/sam/samosa #/net/dials/raid1/sauter/LS49_samosa
export CONDA_ROOT=/net/dials/raid1/mlxd/InDev/STRUMPACK_MPI_DIST
source ${CONDA_ROOT}/miniconda/bin/activate myEnv
mkdir -p ${SAMOSA_ROOT}/${TAG}
cd ${SAMOSA_ROOT}/${TAG}
#data=/net/dials/raid1/sauter/LS49_integ_step4K/int*.pickle \
#data=/net/dials/raid1/sauter/LS49_integ_step4K/int-0-step4K_MPIbatch_0[0-9][0-9][0-9]00.img.pickle \

# Example works, takes 1 hour
# Uses 6554 images, compromise between wall clock time and realistic calculation
# Must set MULTINODE=True to run samosa
# levmar.parameter_flags has not been set. This example treats measurements as fulls & no postrefinement
# To treat spots as partials set either PartialityDeff or PartialityEtaDeff
# Currently refines only scaling factor G by default
# To actually do something, refine any combination of Bfactor Deff Eta Rxy
# But it is burdensome to ask science questions until performance has been improved with Strumpack
#data=/net/dials/raid1/sauter/LS49_integ_step4K/int-0-step4K_MPIbatch_0[0-9][0-9][0-9][0-9][0].img.pickle

export effective_params="d_min=2.3 \
data=/net/dials/raid1/sauter/LS49_integ_step4K/int-0-step4K_MPIbatch_0[0][0-9][0-9][0-9][0-9].img.pickle \
extension=img.pickle \
model=None \
backend=Flex \
scaling.report_ML=True \
pixel_size=0.11 \
nproc=1 \
include_negatives=True \
target_unit_cell=67.2,59.8,47.2,90,110.3,90 \
target_space_group=C2 \
significance_filter.apply=False \
set_average_unit_cell=True \
scaling.algorithm=levmar \
\
levmar.linsolver=cg \
levmar.strumpack.enable=True \
levmar.strumpack.mpi=False \
levmar.strumpack.algorithm=bicgstab \
levmar.strumpack.reordering=scotch \
levmar.strumpack.hss=False \
levmar.strumpack.verbose=True \
\
levmar.compute_cc_half=False \
levmar.sdfac_value=40. \
levmar.termination.step_threshold=0.001 \
levmar.termination.objective_decrease_threshold = 1.E-6 \
memory.shared_array_allocation=40000000 \
scaling.mtz_file=${MERGE_ROOT}/1m2a_fmodel.mtz \
scaling.mtz_column_F=fmodel \
postrefinement.enable=False \
raw_data.errors_from_sample_residuals=True \
raw_data.sdfac_refine=False \
raw_data.error_models.sdfac_refine.random_seed=42 \
min_corr=-1.0 \
output.n_bins=10 \
output.prefix=${trial}"
#levmar.parameter_flags=Bfactor

if [ "${MULTINODE}" == "True" ]; then
  #samosa.join ${effective_params} |& tee samosa_join.jog
  export OMP_NUM_THREADS=16
  (time samosa.scale ${effective_params}) #|& tee samosa_scale.log
  exit
fi
exit
export PHENIX_ROOT=/net/viper/raid1/sauter
#source the phenix build for single node execution:
source ${PHENIX_ROOT}/phenix-dev-2880/build/setpaths.sh

phenix.xtriage ${trial}.mtz scaling.input.xray_data.obs_labels=Iobs |tee xtriage_${trial}.out && \

#xray_data.r_free_flags.file_name=${CONDA_ROOT}/modules/exafel_project/nks/merge118/reflections.mtz \
#xray_data.r_free_flags.label=FreeR_flag \

phenix.refine ${MERGE_ROOT}/1m2a.pdb \
refinement.output.prefix=${trial} \
xray_data.file_name=${trial}.mtz \
xray_data.labels="Iobs" \
refinement.input.xray_data.r_free_flags.file_name=${MERGE_ROOT}/1m2a_flags.mtz \
refinement.input.xray_data.r_free_flags.label=rfree_flags \
main.number_of_macro_cycles=3 optimize_xyz_weight=True \
optimize_adp_weight=True nproc=20 \
ordered_solvent=True ordered_solvent.mode=every_macro_cycle \
refine.adp.individual.anisotropic="element Fe" \
refinement.input.symmetry_safety_check=warning \
refine.strategy="*individual_sites *individual_sites_real_space *rigid_body *individual_adp group_adp tls *occupancies group_anomalous" \
xray_data.force_anomalous_flag_to_be_equal_to=True \
main.wavelength=1.734 --overwrite && \

libtbx.python ${CONDA_ROOT}/modules/exafel_project/nks/map_height_at_atoms.py \
${trial}_001.pdb \
${trial}.mtz \
input.xray_data.labels=Iobs \
input.symmetry_safety_check=warning \
xray_data.r_free_flags.file_name=${trial}_001.mtz \
xray_data.r_free_flags.label=R-free-flags \
selection="element Fe" |tee ${trial}_peakht.log
