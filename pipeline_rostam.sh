#!/bin/zsh

#SBATCH --job-name=vascpp
#SBATCH --output=8N-27_scaling_run02.out
#SBATCH --partition=medusa
#SBATCH --nodes=8
#SBATCH --time=60:00

# set xe

echo -e '\nSCALING RUN 1\n'
srun build/vascpp 27 2

echo -e '\nSCALING RUN 2\n'
srun build/vascpp 27 2

