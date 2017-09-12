#!/bin/bash
for ii in $(ls -d */ | grep r0);
do
        cat $ii/000/out/debug/debug* | grep ",start" | wc -l >> ./total_img.log
        cat $ii/000/out/debug/debug* | grep ",integrate_ok" | wc -l >> ./int_img.log

done
