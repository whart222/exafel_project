from __future__ import absolute_import, print_function, division

message = '''
Goes through all images in a dataset and prints out the information regarding the top N images with the highest RMSD
It also prints out a combined experiment list and reflection table
MPI option available
'''

import os,sys,math
from libtbx.phil import parse
from libtbx.utils import Sorry
from exafel_project.ADSE13_25.command_line.indexing_analytics import params_from_phil
from dxtbx.model.experiment_list import ExperimentList,ExperimentListFactory
from dials.array_family import flex
from libtbx.easy_pickle import load
from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
from scitbx.matrix import col

rmsd_phil_scope = parse('''
  input_path = None
    .type = path
    .help = Path to folder which contains the out folder. For example XXX_rgYYYY in the style of cctbx.xfel
  num_images = 6
    .type = int
    .help = Number of images with highest RMSD whose information to print
  mpi = False
    .type = bool
    .help = Enables usage of mpi part of code
  dump_files = True
    .type = bool
    .help = Dump json and pickle files of the highest RMSD frames.
  dump_basename = high_rmsd
    .type = str
    .help = basename that will be used for dumping experiment list and reflection tables. pickle and json will be \
            appended to the name
''')

def find_rmsd_from_refl_tables(experiments, reflections, num_images):
  rmsd = flex.double()
  all_refl = []
  for ii, expt in enumerate(experiments):
    refl_now = reflections.select(reflections['id'] == ii)
    dR = flex.double()
    for refl in refl_now:
      dR.append((col(refl['xyzcal.mm']) - col(refl['xyzobs.mm.value'])).length())
    rmsd.append(1000.0*math.sqrt(dR.dot(dR)/len(dR)))
  idx_list = list(flex.sort_permutation(rmsd, reverse=True)[0:num_images])
  reqd_expt = ExperimentList()
  reqd_refl = flex.reflection_table()
  for ii,idx in enumerate(idx_list):
    reqd_expt.append(experiments[idx])
    refl = reflections.select(reflections['id'] == idx)
    refl['id'].set_selected(flex.size_t(range(len(refl['id']))),ii)
    reqd_refl.extend(refl)
  return reqd_expt, reqd_refl


def find_rmsd_from_files(filenames, root, num_images, rank=0):
  rmsd = flex.double()
  all_expt = []
  all_refl = []
  for filename in filenames:
    fjson=os.path.join(root, filename)
    experiments = ExperimentListFactory.from_json_file(fjson)
    fpickle = os.path.join(root, filename.split('refined_experiments')[0]+'indexed.pickle')
    reflections = load(fpickle)
    ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
    reflections = ref_predictor(reflections)
    for ii,expt in enumerate(experiments):
      cbf_now = experiments[ii].imageset.get_image_identifier(0).split('/')[-1]
      refl_now = reflections.select(reflections['id'] == ii)
      dR = flex.double()
      for refl in refl_now:
        dR.append((col(refl['xyzcal.mm']) - col(refl['xyzobs.mm.value'])).length())
      rmsd.append(1000.0*math.sqrt(dR.dot(dR)/len(dR)))
      all_expt.append(expt)
      all_refl.append(refl_now)
  # Now sort it
  reqd_expt = ExperimentList()
  reqd_refl = flex.reflection_table()
  idx_list = list(flex.sort_permutation(rmsd, reverse=True)[0:num_images])
  print (idx_list)
  for ii,idx in enumerate(idx_list):
    reqd_expt.append(all_expt[idx])
    refl = all_refl[idx]
    refl['id'].set_selected(flex.size_t(range(len(refl['id']))),ii)
    reqd_refl.extend(refl)
  return (reqd_expt, reqd_refl)


def assign_work(root, mpi=False):
  if mpi:
    if rank == 0:
      iterable = []
      for filename in os.listdir(root):
        if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
        iterable.append(filename)
      print ('done appending from rank 0')
      for i_rank in range(1,size):
        iterable2 = [iterable[i] for i in range(len(iterable)) if (i+i_rank)%size == 0]
        #print ('SENDING STUFF', len(iterable2), iterable2[0], iterable2[-1])
        comm.send(iterable2, dest=i_rank)
      # Rank 0 work
      iterable2 = [iterable[i] for i in range(len(iterable)) if (i+0)%size == 0]
      #print ('RANK 0', len(iterable2), iterable2[0], iterable2[-1])

    else:
      iterable2 = comm.recv(source=0)
  else:
    iterable2 = []
    for filename in os.listdir(root):
      if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
      iterable2.append(filename)
  return iterable2
      #print ('GETTING STUFF', len(iterable2), iterable2[0], iterable2[-1])


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '--h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:], phil_scope=rmsd_phil_scope)
  if params.input_path is None:
    params.input_path='.'
  root = os.path.join(params.input_path, 'out')
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
    iterable2 = assign_work(root, mpi=params.mpi)
    # Make sure num_images to be sorted by each rank is less than len(iterable2)
    if params.num_images <= len(iterable2):
      num_images = params.num_images
    else:
      num_images = len(iterable2)
    results = find_rmsd_from_files(iterable2, root, num_images, rank=rank)
    results = comm.gather(results, root=0)
    comm.barrier()
    if rank == 0:
    # stitch together refl tables and experiment lists
      all_experiments = ExperimentList()
      all_reflections = flex.reflection_table()
      expt_counter = -1
      for result in results:
        expt_list,refl_list = result
        for ii,expt in enumerate(expt_list):
          expt_counter +=1
          all_experiments.append(expt)
          refl = refl_list.select(refl_list['id'] == ii)
          refl['id'].set_selected(flex.size_t(range(len(refl['id']))),expt_counter)
          all_reflections.extend(refl)
      experiments, reflections = find_rmsd_from_refl_tables(all_experiments, all_reflections, num_images)
      if params.dump_files:
        from dxtbx.model.experiment_list import ExperimentListDumper
        from libtbx.easy_pickle import dump
        el_dumper = ExperimentListDumper(experiments)
        el_dumper.as_json('high_rmsd.json')
        dump('high_rmsd.pickle', reflections)



  else:
    iterable2 = assign_work(root, mpi=params.mpi)
    if params.num_images <= len(iterable2):
      num_images = params.num_images
    else:
      num_images = len(iterable2)
    results = find_rmsd_from_files(iterable2, root, num_images, rank=0)
    #from IPython import embed; embed(); exit()
    if True:
    # stitch together refl tables and experiment lists
      all_experiments = ExperimentList()
      all_reflections = flex.reflection_table()
      expt_list,refl_list = results
      for ii,expt in enumerate(expt_list):
        all_experiments.append(expt)
        refl = refl_list.select(refl_list['id'] == ii)
        refl['id'].set_selected(flex.size_t(range(len(refl['id']))),ii)
        all_reflections.extend(refl)
      experiments, reflections = find_rmsd_from_refl_tables(all_experiments, all_reflections, num_images)
      if params.dump_files:
        from dxtbx.model.experiment_list import ExperimentListDumper
        from libtbx.easy_pickle import dump
        el_dumper = ExperimentListDumper(experiments)
        el_dumper.as_json(params.dump_basename +'.json')
        dump(params.dump_basename+'.pickle', reflections)
