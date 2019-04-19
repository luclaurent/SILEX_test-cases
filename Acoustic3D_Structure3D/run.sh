export OMP_NUM_THREADS=1
mpirun -np 20 python parametricStudyXY.py >> log.out 2>&1 &
