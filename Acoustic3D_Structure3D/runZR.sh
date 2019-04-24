export OMP_NUM_THREADS=1
mpirun --hostfile hostfile  -np 16  python parametricStudyZR.py >> logZR.out 2>&1 &
