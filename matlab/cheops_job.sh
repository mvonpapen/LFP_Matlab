#!/bin/bash -l

#SBATCH --job-name myprog
#SBATCH --cpus-per-tasks=4
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --account=AG-Saur

module load matlab

myprog="/scratch/matlab/Granger/timeavg.m"

# start matlab from shell
time matlab -nodisplay -nodesktop -nosplash -nojvm \
	    -r "run $myprog"
