#!/bin/bash

#PBS -o /home/aeurai/jobs
#PBS -e /home/aeurai/jobs
#PBS -lnodes=1 -lwalltime=00:05:00

# load necessary modules
module load matlab

# run the plotting script
matlab -nodesktop -nodisplay -r "disp(pwd); plot_all; end; exit"
