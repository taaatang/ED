#!/bin/bash -l

#SBATCH --job-name=k6by6_kIndex
#SBATCH -C knl
#SBATCH -L SCRATCH,project 
#SBATCH -q regular
#SBATCH --time=48:00:00
#SBATCH --nodes=50
#SBATCH --ntasks=200
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=300
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tatang@stanford.edu
export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
#export OMP_DISPLAY_AFFINITY=TRUE
#export OMP_AFFINITY_FORMAT="host=%H, pid=%P, thread_num=%n, thread affinity=%A"
srun --cpu-bind=cores ./kGroundState