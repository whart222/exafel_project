EXP=$(cat webUI.json | python -c 'import json,sys;obj=json.load(sys.stdin);print obj["experiment_name"];’)
RUNNO=$(cat webUI.json | python -c 'import json,sys;obj=json.load(sys.stdin);print obj["run_num"];')
RUNFOLDER=$(printf “r%04d” $RUNNO)
TRIALNO=0 # Fireworks should keep track of what TRIALNO to use
TRIALFOLDER=$(printf “%03d” $TRIALNO)

cctbx.xfel.xtc_process \
input.experiment=$EXP \
input.run_num=$RUNNO \
output.logging_dir=dials/$RUNFOLDER/$TRIALFOLDER/stdout \
output.output_dir=dials/$RUNFOLDER/$TRIALFOLDER/out \
target.phil
