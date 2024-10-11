from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank  = comm.Get_rank()

# To run it:
# mpirun -np 2  python exemple_mpi.py

print 'hello, je suis le numero : ',rank

a1=3.0
b1=6.0

a2=7.0
b2=9.0


if rank==0:
    x1=b1/a1

if rank==1:
    x2=b2/a2
    comm.send(x2, dest=0, tag=11)

if rank==0:
    data = comm.recv(source=1, tag=11)
    print 'x1=',x1
    print 'x2=',data
