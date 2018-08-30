#!/bin/bash
#Extract the anom peak height from the resulting files

for ii in $(find $PWD -type f -iname "tol*ErrRes*avg*.log" | grep -v mark0 | grep -v graphs | grep -v peak | grep -v 001);
do
  echo $(basename $ii);
  sed '9q;d' $ii;
done
for ii in $(find $PWD -type f -iname "*peak*log");
do
  echo $(basename $ii);
  tail -n 1 $ii
done
