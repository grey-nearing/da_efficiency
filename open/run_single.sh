#!/bin/bash

# screen report
echo "running site: $1" 

# go into directory
cd "./runtime_dirs/site_$1"

# run experiment
./main.exe 

# gather mean results
cp open_mean.out ../../site_data/open_$1.out 

# delete runtime directory to save space
cd ../..
/bin/rm -r ./runtime_dirs/site_$1

# screen report
echo "finished site: $1" 


 


