#!/bin/bash

for (( o = 1 ; o <= 11 ; o++ )); do
 for (( m = 1 ; m <= 15 ; m++ )); do
#  factor=$(echo $f | awk '{printf "%4.3f\n",$1*0.1}')
  sh run_single_obs_cov.sh $1 $o $m
 done
done

/bin/rm -r ./runtime_dirs/site_$1

