Installation script:

https://github.com/ExaFEL/adse13_161/blob/master/summit/dials_py36/install_py36_dials.sh

After the installation please overwrite the modules/ccp4io/libccp4/ccp4/library_file.c file with the file in this git folder. This is necessary to prevent ccp4 code from swapping bytes in mtz files on Summit. Then re-run make.

Submitting a merging batch job:

bsub merge_LD91.lsf

For super-large datasets like PSII, first run calculate_file_load.sh to determine the minimum required number of nodes.
Then in your merging batch script request at least that many nodes. Also, in your merging script:
1. Make sure the input.parallel_file_load parameters are the same as in calculate_file_load.sh
2. If using more than the minimum required number of nodes, be sure to also specify input.parallel_file_load.balance=global to spread the input data across CPU cores
3. Should your load balancing fail with an mpi4py alltoall error, do alltoall iteratively by specifying  input.parallel_file_load.balance_mpi_alltoall_slices=2 or larger 
