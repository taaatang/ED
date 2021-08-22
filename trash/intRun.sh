export OMP_NUM_THREADS=16
srun -N 2 -n 8 -c 16 --mem-per-cpu 1000 -o kGroundState7by4N2.out ./kGroundState 7 4
