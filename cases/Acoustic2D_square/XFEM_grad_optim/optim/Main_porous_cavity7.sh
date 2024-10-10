#! /bin/sh
#!/usr/bin/python
nbProc=$3
angle=$1
positionX=$2
# Make the mesh for the given parameters
python3.4 Main_make_mesh_porous_cavity7.py $angle $positionX

# launch the computation
export OPENBLAS_NUM_THREADS=1
# here use 20 proc. to compute 20 freq. steps in parallel
mpirun -np $nbProc python3.4 Main_classic_fluid_porous7.py 

# plot the frf
#python3.4 Plot_frf_fluid_porous7.py

# exemple of use:
#./Main_porous_cavity7.sh 10.0 2.0

