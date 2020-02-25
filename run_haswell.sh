#!/bin/bash -l

#SBATCH --job-name=4x7spec
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -L SCRATCH,project 
#SBATCH -t 24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=32
# #SBATCH --mem-per-cpu=1500
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tatang@stanford.edu
export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
# export OMP_DISPLAY_AFFINITY=TRUE
# export OMP_AFFINITY_FORMAT="host=%H, pid=%P, thread_num=%n, thread affinity=%A"
srun --cpu-bind=cores ./Spectra.out