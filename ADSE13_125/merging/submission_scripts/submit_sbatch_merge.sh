merge_script="#!/bin/bash -l
\n
\n#SBATCH -C knl
\n#SBATCH -N 1
\n#SBATCH -q regular
\n#SBATCH -J m_$1_$2
\n#SBATCH -t 24:00:00 
\n#SBATCH -A m2859 
\n
\nsource /global/u2/a/asmit/cctbx.xfel/build/setpaths.sh
\n
\nfor i in `seq -f "%04g" $1 $2`; do ./merge.sh \$i; done"


#cat <<$script>> out
echo -e $merge_script > sbatch_merge.sh
sbatch sbatch_merge.sh
sleep 5
