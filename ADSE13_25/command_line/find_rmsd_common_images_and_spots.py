from __future__ import absolute_import,print_function, division

message = ''' Script to calculate RMSDs resulting from different modes of indexing (or even just one mode). If multiple folders are specified the common set of images+observed spots are taken and used for RMSD calculation. If single folder is specified, all the images are used
'''

import sys, os, math
from dials.array_family import flex
from libtbx.utils import Sorry # implicit import
from exafel_project.ADSE13_25.command_line.indexing_analytics import params_from_phil
from libtbx.phil import parse

rmsd_phil_scope = parse('''
  input_path = None
    .multiple = True
    .type = path
    .help = input_path should be like that used at LCLS i.e full path to run_number/trial_rg \
            Assumes folder structure \
            Input Path    -----> out -----> debug folder \
                          -----> pickle and json files \
                          -----> stdout \
            You can specificy multiple paths. In that case, the unit cells and RMSD info will be \
            only for the common set of images indexed
  mpi = False
    .type = bool
    .help = If True, mpi can be used for running the program
''')

def get_common_set(roots):
  ''' Function to get common set of images+spots indexed in multiple folders. Based on CBF filenames
      ts_from_cbf if True is much faster than reading from json files
      Returns a dictionary of timestamps with common spots recorded for each ts'''
  from dxtbx.model.experiment_list import ExperimentListFactory
  from libtbx.easy_pickle import load
  cbf = {}
  for root in roots:
    cbf[root] = []
    if True:
      for filename in os.listdir(root):
        if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
        explist=ExperimentListFactory.from_json_file(os.path.join(root,filename), check_format=False)
        for exp in explist:
          cbf[root].append(exp.imageset.get_image_identifier(0).split('/')[-1])
  # Now take intersection
  common_set = set(cbf[roots[0]])
  for ii in range(1, len(roots)):
    common_set = common_set.intersection(set(cbf[roots[ii]]))
  common_set=list(common_set)

  # Now that we have common set of images, let's get the common set of spots (observed)
  obs_spots = {}
  for root in roots:
    for filename in os.listdir(root):
      if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
      explist=ExperimentListFactory.from_json_file(os.path.join(root,filename), check_format=False)
      fpickle = os.path.join(root, filename.split('refined_experiments')[0]+'indexed.pickle')
      reflections = load(fpickle)
      for ii,exp in enumerate(explist):
        ts=exp.imageset.get_image_identifier(0).split('/')[-1]
        if ts not in common_set: continue
        if ts not in obs_spots:
          obs_spots[ts] = {}
        refl=reflections.select(reflections['id']==ii)
        obs_spots[ts][root] = list(refl['xyzobs.mm.value'])

  common_image_and_spot_set = {}
  for ts in common_set:
    common_image_and_spot_set[ts] = set(obs_spots[ts][roots[0]])
    for ii in range(1, len(roots)):
      common_image_and_spot_set[ts]=common_image_and_spot_set[ts].intersection(set(obs_spots[ts][roots[ii]]))
    common_image_and_spot_set[ts] = list(common_image_and_spot_set[ts])
  return common_image_and_spot_set

def get_rmsd_stats(filenames,root, rank=0,common_set=None):
  ''' Provide a folder with refined_experiment.json and indexed pickle files. It will calculate the RMSD'''

  from dxtbx.model.experiment_list import ExperimentListFactory
  from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
  from libtbx.easy_pickle import load
  from scitbx.matrix import col

  dR = flex.double()
  if common_set is not None:
    print ('Using %d common set images to report RMSD statistics'%(len(common_set)))
  for filename in filenames:
    fjson=os.path.join(root, filename)
    experiments = ExperimentListFactory.from_json_file(fjson, check_format=False)
    fpickle = os.path.join(root, filename.split('refined_experiments')[0]+'indexed.pickle')
    reflections = load(fpickle)
    ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
    reflections = ref_predictor(reflections)
    for ii,expt in enumerate(experiments):
      if common_set is not None:
        cbf_now = experiments[ii].imageset.get_image_identifier(0).split('/')[-1]
        if cbf_now not in common_set: continue
      refl=reflections.select(reflections['id']==ii)
      for entry in refl:
        if common_set is not None:
          if entry['xyzobs.mm.value'] not in common_set[cbf_now]: continue
        dR.append((col(entry['xyzcal.mm']) - col(entry['xyzobs.mm.value'])).length())
  return dR

def run(params, root, common_set=None):
  iterable2 = []
  for filename in os.listdir(root):
    if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
    iterable2.append(filename)
  if params.mpi:
    try:
      from mpi4py import MPI
    except ImportError:
      raise Sorry("MPI not found")
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print (rank, size)
    # get rmsd info
    iterable2 = [iterable2[i] for i in range(len(iterable2)) if (i+rank)%size == 0]
    results2 = get_rmsd_stats(iterable2, root,rank=rank,common_set=common_set)
    results2 = comm.gather(results2, root=0)
    if rank != 0: return
  else:
    results2 = [get_rmsd_stats(iterable2, root, rank=0, common_set=common_set)]

  dR = flex.double()
  info_list = []
  info = []
  for ii,r in enumerate(results2):
    dR.extend(r)
  print ('Total RMSD i.e calc - obs for Bragg spots (um) = ', 1000.0*math.sqrt(dR.dot(dR)/len(dR)))

if __name__=='__main__':
  params = params_from_phil(sys.argv[1:], phil_scope=rmsd_phil_scope)
  if len(params.input_path) == 0:
    params.input_path=['.']
  roots = []
  for input_path in params.input_path:
    roots.append(os.path.join(input_path, 'out'))
  common_set=None
  if len(roots) > 1:
    common_set = get_common_set(roots)
  for root in roots:
    run(params,root, common_set=common_set)
