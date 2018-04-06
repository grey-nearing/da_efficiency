#!/bin/bash

#for (( f = 1 ; f <= 20 ; f++ )); do
 for (( s = 1 ; s <= 129 ; s++ )); do

#  factor=$(echo $f | awk '{printf "%4.3f\n",$1*0.1}')
  echo sh run_site_obs_cov.sh $s # $factor

# done
done > 'execfile.slurm'


