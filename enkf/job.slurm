#!/bin/bash -x

#SBATCH -N 1
#SBATCH -n 40
#SBATCH --job-name=e_xx_yy
#SBATCH --output=report_xx_yy.slurm
#SBATCH --time=36:00:00
#SBATCH -q defq
#SBATCH --exclude=gpu01

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

echo xx yy
sh execfile.sh xx yy
wait

# done
exit 0


