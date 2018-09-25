for i in `seq 175 255`
do
  for j in `ls ../r0${i}/001_rg037/out/*integrated_experiments.json`
  do
    #echo $j
    run=_r0${i}
    prefix=../r0${i}/001_rg037/out/
    echo $prefix
    flink=${j#"$prefix"}
    echo $flink
    p=8
    flink2="${flink:0:p}${run}${flink:p}"
    echo $flink2
    ln -s $j $flink2
  done

done

for i in `seq 175 255`
do
  for j in `ls ../r0${i}/001_rg037/out/*integrated.pickle`
  do
    #echo $j
    run=_r0${i}
    prefix=../r0${i}/001_rg037/out/
    echo $prefix
    flink=${j#"$prefix"}
    echo $flink
    p=8
    flink2="${flink:0:p}${run}${flink:p}"
    echo $flink2
    ln -s $j $flink2
  done

done
