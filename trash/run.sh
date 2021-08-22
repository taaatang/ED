#!/bin/bash -l

#SBATCH -q regular
#SBATCH -N 2
#SBATCH --tasks-per-node 4
#SBATCH -J k7by4
#SBATCH -L SCRATCH,project 
#SBATCH -C haswell
#SBATCH --mem-per-cpu 1000
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tatang@stanford.edu
export OMP_NUM_THREADS=8
srun -n 8 -c 16 -o kGroundState7by4.out ./kGroundState 7 4


