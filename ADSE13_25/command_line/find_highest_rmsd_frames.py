from __future__ import absolute_import, print_function, division

message = '''
Goes through all images in a dataset and prints out the information regarding the top N images with the highest RMSD
'''

import os,sys
from libtbx.phil import parse
from libtbx.utils import Sorry
from exafel_project.ADSE13_25.command_line.indexing_analytics import params_from_phil

rmsd_phil_scope = parse('''
  input_path = None
    .type = path
    .help = Path to folder which contains the out folder. For example XXX_rgYYYY in the style of cctbx.xfel
  num_images = 100
    .type = int
    .help = Number of images with highest RMSD whose information to print
  mpi = False
    .type = bool
    .help = Enables usage of mpi part of code
''')


def find_rmsd(params):
  return

def assign_work(params, root, mpi=False)
  if mpi:
    if rank == 0:
      iterable = []
      for filename in os.listdir(root):
        if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
        iterable.append(filename)
      print ('done appending from rank 0')
      for i_rank in range(1,size):
        iterable2 = [iterable[i] for i in range(len(iterable)) if (i+i_rank)%size == 0]
        comm.send(iterable2, dest=i_rank)
      
    else:
      iterable2 = comm.recv(source=0)


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '--h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:], phil_scope=venn_phil_scope)
  if params.input_path is None:
    params.input_path='.'
  root = os.path.join(input_path, 'out')
  if params.mpi:
    try:
      from mpi4py import MPI
    except ImportError:
      raise Sorry("MPI not found")
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print (rank, size)
    # Assign work to ranks
    if rank == 0:
      iterable = []
      for filename in os.listdir(root):
        if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
        iterable.append(filename)
      print ('done appending from rank 0')
    else:
      iterable = None
    iterable = comm.bcast(iterable, root=0)
    iterable = [iterable[i] for i in range(len(iterable)) if (i+rank)%size == 0]
    # Now get RMSD info
    results = find_rmsd_stats(iterable, root,rank=rank)
    results = comm.gather(results, root=0)
    if rank != 0: return
  else:
    results = [get_hits_and_indexing_stats(iterable, debug_root)]
    results2 = [get_uc_and_rmsd_stats(iterable2, root, rank=0, common_set=common_set)]

