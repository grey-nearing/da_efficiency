#!/bin/bash

da_flag=-50

for (( i = 1 ; i <= 129; i++ )); do
   echo $i

   rm -rf runtime_dirs/site_$i
   mkdir runtime_dirs/site_$i

   cp ./setup_dir/init.txt runtime_dirs/site_$i/.
   echo $da_flag > runtime_dirs/site_$i/da_flag.txt
   cp ./setup_dir/main.exe runtime_dirs/site_$i/.
   cp ./setup_dir/num_times.txt runtime_dirs/site_$i/.
   head -$i ../data/scan/SCAN_Site_List.txt | tail -1 > DelThisSCAN_Site.txt
   cut -f 2 DelThisSCAN_Site.txt > lat.txt 
   cut -f 3 DelThisSCAN_Site.txt > lon.txt 
   offset=$(cut -f 4 DelThisSCAN_Site.txt)
   cat lat.txt lon.txt >  runtime_dirs/site_$i/lat_lon.txt
   rm -rf lat.txt lon.txt DelThisSCAN_Site.txt

   for (( j = 1 ; j <= -$da_flag; j++ )); do
     ln -s ../../../data/forcing/perturbed/forcing_${i}_${j}.txt ./runtime_dirs/site_$i/forcing_$j.txt
   done

   cp ../data/parms/extract_parms/site_data/parms_$i.txt ./runtime_dirs/site_$i/parms.txt
   cp ../data/parms/extract_parms/site_data/time_parms_$i.txt ./runtime_dirs/site_$i/time_parms.txt
   cp ../data/lprm/site_data/lprm_cdf_$i.txt ./runtime_dirs/site_$i/obs.txt
   cp ../data/lprm/site_data/cdf_sig_$i.txt ./runtime_dirs/site_$i/obs_cov.txt
   cp ./setup_dir/plant_init.txt runtime_dirs/site_$i/.
   cp ./setup_dir/soil_init.txt runtime_dirs/site_$i/. 
   cp ./setup_dir/startdate.txt runtime_dirs/site_$i/.

#   let startdate=200012311600+$((8-$offset))\*100
#   echo $startdate > runtime_dirs/site_$i/startdate.txt

   echo "finished site: $i" 

done > report.setup_dirs   

 


