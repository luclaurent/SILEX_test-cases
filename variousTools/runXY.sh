export OMP_NUM_THREADS=1
mpirun -np 20 python parametricStudyXY.py >> logXY.out 2>&1 &
