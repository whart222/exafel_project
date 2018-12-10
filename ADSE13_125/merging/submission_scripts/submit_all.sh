./submit_sbatch_merge.sh 33 1000
# 1000 - 34000
step=1000
for i in `seq 1 33`
do
  count=$(($i * $step))
  run1=$(($count+1))
  j=$(($i+1))
  run2=$(($j * $step))
  echo $run1 $run2
  ./submit_sbatch_merge.sh $run1 $run2
done 
./submit_sbatch_merge.sh 34001 34816
