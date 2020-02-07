Installation script:

https://github.com/ExaFEL/adse13_161/blob/master/summit/dials_py36/install_py36_dials.sh

Submitting a merging batch job:

bsub merge_LD91.lsf

For super-large datasets like PSII, first run calculate_file_load.sh to determine the minimum required number of nodes.
Then in your merging batch script:
1. Request at least that many nodes
2. Make sure that your input.parallel_file_load parameters are the same as in calculate_file_load.sh
3. If requesting more than the minimum required number of nodes, be sure to also specify input.parallel_file_load.balance to spread the input data across CPU cores
4. If your load balancing fails with an mpi4py alltoall error, you will have to do alltoall iteratively by specifying: input.parallel_file_load.balance_mpi_alltoall_slices 
