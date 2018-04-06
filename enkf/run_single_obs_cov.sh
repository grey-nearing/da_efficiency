#!/bin/bash

# screen report
echo "running site: $1 -- obs cov factor = $2 -- state pert factor = $3" 

# read in obs_pert and state_pert from master list file
obsPert=$(head -"$2" 'obs_perts_master.txt' | tail -1)
statePert=$(head -"$3" 'state_perts_master.txt' | tail -1)

# go into directory
cd "./runtime_dirs/site_$1"

# change obs_cov file
#cp ../../../data/lprm/site_data/cdf_sig_$1.txt obs_cov.txt
#typeset -r obsCov=$(cat 'obs_cov.txt')
#typeset -r factor=$obsPert
#obsCovFac=$(echo $obsCov $factor | awk '{printf "%7.5f\n",$1*$2}')
#echo $obsCovFac > 'obs_cov.txt'
echo $obsPert > 'obs_cov.txt'

# set state pert for this run
echo $statePert $statePert $statePert $statePert > 'state_pert.txt'

# run experiment
./main.exe #> Screenout

# gather mean results
mv enks_mean.out ../../site_data/enkf_$1_$2_$3.out 
mv back_mean.out ../../site_data/back_$1_$2_$3.out 

# gatehr ensemble results
# for (( j = 1 ; j <= $da_flag; j++ )); do
#  cp back_$j.out ../site_data/back_${i}_${j}.out 
#  cp enks_$j.out ../site_data/enks_${i}_${j}.out 
# done

# screen report
echo "finished site: $1 -- obs cov factor = $2 -- state pert factor = $3" 


 


