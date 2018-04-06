#!/bin/bash

#Nsites=129; 383
#SY: NOTE: change range of ss in ss loop below as required

for ss in {1..383}
do
 echo $ss
 cp ../gather_parms/site_data/sparse/parms_${ss}.txt ./parms.txt
 ./driver.exe namelist
 cp parms.out ../site_data/sparse/parms_${ss}.txt
 cp time_parms.txt ../site_data/sparse/time_parms_${ss}.txt
done
\rm parms.txt parms.out time_parms.txt



