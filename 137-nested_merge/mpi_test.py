from __future__ import division
from mpi4py import MPI

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

print "Rank %d of %d [%s]."%(rank, size, name)
