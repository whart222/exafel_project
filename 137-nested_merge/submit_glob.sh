#!/bin/bash
shopt -s extglob
export TARDIR=$1 #Pass in the directory with the tar files:   eg  /my/scratch/dir/TAR_95-114
export OUTDIR=$2 #Pass in the directory for burst-buffer data to be staged out:   /my/scratch/dir/bb_out
RUNARRAY=() #Array to hold the structure of the files to be processed

#Iterate over the given runs and generate the TAR globs
echo "Generating expanded globs for TAR file paths"
for ii in $(seq -w 95 114); do
  RUNARRAY+=($ii)
  str=""
  for jj in ${RUNARRAY[@]}; do
    str+=$jj,
  done
  str=$(echo $str | sed s'/.$//')
  if [ "${#RUNARRAY[@]}" == "1" ]; then
    g=$(eval "echo \"r0${str}*.tar\" ")
  else
    g=$(eval "echo \"r0{${str}}*.tar\" ")
  fi
  echo $g
  
  #Create a submission script with the appropriate arguments
  cat sbatch_merge.sh | sed "s~<tag_template>~r0${ii}~g" | sed "s~<glob_template>~${g}~g" | sed "s~<tar_template>~${TARDIR}~g"  | sed "s~<data_out_template>~${OUTDIR}~g" > sbatch_merge_r0${ii}.sh

  #Submit the newly created script to the queue
 # sbatch sbatch_merge_r0${ii}.sh $PWD
done
