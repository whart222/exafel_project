from __future__ import absolute_import, print_function, division

message = ''' Script to plot residual vectors given pairs of experiment json and reflection pickle files.
Uses certain components of dials_refinement_preceding_integration.py in cctbx_project/rstbx'''

import os,sys,math
from libtbx.phil import parse
from libtbx.utils import Sorry
from exafel_project.ADSE13_25.command_line.indexing_analytics import params_from_phil
from dxtbx.model.experiment_list import ExperimentList,ExperimentListFactory
from dials.array_family import flex
from libtbx.easy_pickle import load
from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
from scitbx.matrix import col

residual_vector_plot_phil_scope = parse('''
  input_path = None
    .type = path
    .help = Path to where reflection tables and json files are located
  basename = None
    .type = str
    .help = If you want to plot residual vectors of a certain json-pickle combo in a directory\
            full of those files, specify the basename i.e basename.json and basename.pickle
  show_plot = True
    .type = bool
    .help = Flag to indicate if plots should be displayed on screen
  show_residual_scatter_plot = False
    .type = bool
    .help = Flag to indicate if residual scatter plot is to be displayed
  show_residual_map_plot = True
    .type = bool
    .help = Flag to indicate if residual map plot is to be displayed. 
''')




# Read in an experiment list and reflection table
def plot_residual_vectors(experiment_list_files, reflection_table_files, params):
  correction_vectors_provisional = []
  indexed_pairs_provisional = []
  for fjson, fpickle in zip(experiment_list_files, reflection_table_files):
    experiments = ExperimentListFactory.from_json_file(fjson, check_format=False)
    reflections = load(fpickle)
    ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
    reflections = ref_predictor(reflections)
    for ii,expt in enumerate(experiments):
      #cbf_now = experiments[ii].imageset.get_image_identifier(0).split('/')[-1]
      FWMOSAICITY=2*expt.crystal.get_half_mosaicity_deg()
      DOMAIN_SZ_ANG=expt.crystal.get_domain_size_ang()
      refl_now = reflections.select(reflections['id'] == ii)
      this_setting_matched_indices = refl_now["miller_index"]
      hkllist = refl_now['miller_index']
      for j,item in enumerate(this_setting_matched_indices):
        # Every miller index has a prediction. So assuming this_setting_index same as j 
        this_setting_index = hkllist.first_index(item)
        Match = dict(spot=j,pred=this_setting_index)
        indexed_pairs_provisional.append(Match)
      for refl in refl_now:
        vector = col((refl['xyzobs.px.value'][0]-refl['xyzcal.px'][0], refl['xyzobs.px.value'][1]-refl['xyzcal.px'][1]))
        correction_vectors_provisional.append(vector)

  if params.show_plot and params.show_residual_scatter_plot:
    from matplotlib import pyplot as plt
    fig = plt.figure()
    for cv in correction_vectors_provisional:
      plt.plot([cv[1]],[-cv[0]],"r.")
    plt.axes().set_aspect("equal")


  PX = reflections["xyzobs.px.value"]
  if params.show_plot and params.show_residual_map_plot:
    from matplotlib import pyplot as plt
    fig = plt.figure()
    for match,cv in zip(indexed_pairs_provisional,correction_vectors_provisional):
      plt.plot([PX[match["spot"]][1]],[-PX[match["spot"]][0]],"r.")
      plt.plot([reflections['xyzcal.px'][match["pred"]][1]],[-reflections['xyzcal.px'][match["pred"]][0]],"g.")
      plt.plot([PX[match["spot"]][1], PX[match["spot"]][1] + 10.*cv[1]],
             [-PX[match["spot"]][0], -PX[match["spot"]][0] - 10.*cv[0]],'r-')
    if True:
      from rstbx.apps.stills.util import residual_map_special_deltapsi_add_on
      residual_map_special_deltapsi_add_on(
        reflections = reflections,
        matches = indexed_pairs_provisional, experiments=experiments,
        hkllist = hkllist,
        predicted = reflections['xyzcal.mm'], plot=plt, eta_deg=FWMOSAICITY, deff=DOMAIN_SZ_ANG
        )
#  plt.xlim([0,detector[0].get_image_size()[1]])
#  plt.ylim([-detector[0].get_image_size()[0],0])
#  plt.title(" %d matches, r.m.s.d. %5.2f pixels"%(len(correction_vectors_provisional),math.sqrt(flex.mean(c_v_p_flex.dot(c_v_p_flex)))))
    plt.axes().set_aspect("equal")
  if params.show_plot and (params.show_residual_scatter_plot or params.show_residual_map_plot):
    plt.show()

if __name__ == '__main__':
  params = params_from_phil(sys.argv[1:], phil_scope=residual_vector_plot_phil_scope) 
  if params.input_path is None:
    root=os.getcwd()
  else:
    root=params.input_path
  experiment_list_files = []
  reflection_table_files = []
  for filename in os.listdir(root) :
    if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
    experiment_list_files.append(os.path.join(root, filename))
    reflection_table_files.append(os.path.join(root, filename.split('refined_experiments')[0]+'indexed.pickle'))
  plot_residual_vectors(experiment_list_files, reflection_table_files, params)
