#!/bin/bash

sed "s/xx/$1/g" job.slurm > temp_job.slurm
sed -i "s/yy/$2/g" temp_job.slurm
sbatch temp_job.slurm
rm temp_job.slurm

