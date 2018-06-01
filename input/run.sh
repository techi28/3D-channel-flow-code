#!/bin/bash
#SBATCH --job-name=channel
#SBATCH --mem=4G
#SBATCH --nodes=2
#SBATCH --partition=cn
#SBATCH --ntasks-per-node=8
#SBATCH --output=channel-%J.out
#SBATCH --time=70:00:00

echo "Starting..."
echo "Environment:"
env | sort > env.log
mpirun -np 16 ./nsmp3D > log.txt
echo "Stopping..."
