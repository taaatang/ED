#!/bin/bash
#SBATCH --job-name=2x2eq
#SBATCH --partition=iric,normal
#SBATCH --time=0-01:00:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
##SBATCH --mem-per-cpu=4000
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=taaatang@gmail.com
#export OMP_NUM_THREADS=4
mpirun -n 1 /home/users/tatang/project/PD/time/ED/build/main.out