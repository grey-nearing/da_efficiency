#!/bin/bash -x

#SBATCH -N 1
#SBATCH -n 40
#SBATCH --job-name=open
#SBATCH --output=report.slurm
#SBATCH --time=12:00:00
#SBATCH -q defq

# memory
ulimit -s unlimited

## modules
#source /usr/share/modules/init/bash
#module purge
#module load comp/intel-16.0.3.210
#
## linking
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/discover/nobackup/projects/lis/libs/netcdf/4.3.3.1_intel-14.0.3.174_sp3/lib
#envVar="{'LD_LIBRARY_PATH': '$LD_LIBRARY_PATH'}"
#
## location
#cd /discover/nobackup/gnearing/projects/da_efficiency/open
#
## pods
#/usr/local/other/PoDS/PoDS/pods.py -e " $envVar " -x /discover/nobackup/gnearing/projects/da_efficiency/open/execfile.slurm -n 10

sh execfile.slurm
wait

# done
exit 0


