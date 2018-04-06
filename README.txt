Welcome to the REAMDE file!  

How To Run:

1) Compile Noah-MP model:
	1a) Need a ifort compiler. On NCCS Discover, this can be accessed by sourcing the <environment_setup.sh> file. Note that this purges all modules and changes $PATH.
	1b) In the <model_code> subdirectory run 'make clean' and 'make'.

2) Run the open-loop simulations:
	2a) Naviage to the <open> subdirectory.
	2b) Create runtime subdirectories by running <setup_dirs.sh> script.
	2c) Edit the hard-paths and header information in the <job.slurm> SLURM/PODS script [*].
	2d) Run the open-loop simulations by submitting the <job.script> SLURM file.

3) CDF-match the LPRM observations:
	3a) Navigate to the <data/lprm> subdirectory.
	3b) Run the <cdf_matching.m> matlab script.

4) Run the EnKF-loop simulations:
	4a) Navigate to the <enkf> subdirectory.
	4b) Create runtime subdirectories by running <setup_dirs.sh> script.
	4c) Edit the hard-paths and header information in the <job.slurm> SLURM/PODS script [*].
	4d) Run the open-loop simulations by submitting the <job.script> SLURM file.

5) Run the various analysis MatLab scripts in the <analysis> subdirectory to create teh various figures in the paper.
	5a) <make_scan_map> makes Figure 3 - showing locations of SCAN evaluation data stations.
	5b) <main_convergence> makes Figure 4 - showing convergence of the IT metrics as a function of sample size.
	5c) <main_efficiency> makes Figure 5 & Table 3 - showing main efficiency metrics.
	5d) <main_decomposition> makes Figure 6 & Table 4 - showing decomposition metrics. 
	5e) <sig_vs_info> makes Figure 7 - comparing Gaussian vs. Information-Theory statistics.

[*] On computers other than NCCS Discover, the user will likely have to set up their own job scheduling routines. This one is specifically designed for SLURM with PODS, and uses modules availalbe on Discover.


