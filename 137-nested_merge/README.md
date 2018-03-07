To generate sbatch submission scripts for the nested merge of LD91 runs [[[95],96],..,114] use:
```bash
./submit_glob.sh <path_to_tars> <output_dir_for_data>
```

The script will generate sbatch scripts for each merge data set, named as `sbatch_merge_r0${ii}.sh` where `${ii}` is the index of the highest run to merge.
The script will automatically generate and submit the merges to the queue.


A sample run-script, for Tar files located @ `TARDIR=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/TAR_95-114/` and output directory for data @ `OUTDIR=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/bb/`


```bash
export TARDIR=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/TAR_95-114/
export OUTDIR=/global/cscratch1/sd/mlxd/sept_sprint/merge_multi/bb/
./submit_glob.sh $TARDIR $OUTDIR
```
