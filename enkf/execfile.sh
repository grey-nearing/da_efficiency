#!/bin/bash

nsites=$(($2-$1))
count=0
for (( site = $1 ; site <= $2 ; site++ )); do
  sh run_site_obs_cov.sh $site & 
  count=$(($count+1))
  concurrent=$(($count % 30))
  if [ $count -eq 29 ]
  then
    wait
  fi
done
wait
