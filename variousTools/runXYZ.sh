export OMP_NUM_THREADS=1
mpirun -np 20 python parametricStudyXYZ.py >> logXYZ.out 2>&1 &
