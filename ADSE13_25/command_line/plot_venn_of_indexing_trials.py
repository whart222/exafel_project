from __future__ import absolute_import, division, print_function

message = '''
Plots venn diagram of indexing results from multiple folders
This module requires a module called matplotlib_venn.
Please do pip install matplotlib_venn if you don't have it in your python installation
'''

from matplotlib import pyplot as plt
import os,sys
from libtbx.phil import parse
from libtbx.utils import Sorry
from exafel_project.ADSE13_25.command_line.indexing_analytics import params_from_phil

venn_phil_scope = parse('''
  input_path = None
    .multiple = True
    .type = path
    .help = Path to folders whose indexing results are to be compared and plotted as Venn diagram \
            Using cctbx.xfel gui convention, you can provide the path to XXX_rgYYYY \
            Assumes following directory structure \
            XXX_rgYYYY      \
              -----> out     \
            Path can be relative to where the script is being run from or absolute paths
  show_plot = True
    .type = bool
    .help = Flag to indicate if plot should be displayed
  pickle_plot = False
    .type = bool
    .help = If True, will pickle matplotlib session so that it can be opened later for analysis/viewing \
            https://stackoverflow.com/questions/29160177/matplotlib-save-file-to-be-reedited-later-in-ipython
  pickle_filename = venn_plot.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
''')


def get_indexed_ts(roots):
  ''' Function to get timestamps of images indexed in provided folders. Based on CBF image identifier '''
  from dxtbx.model.experiment_list import ExperimentListFactory
  cbf = []
  for root in roots:
    idx_filelist = []
    for filename in os.listdir(root):
      if 'refined_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue
      explist=ExperimentListFactory.from_json_file(os.path.join(root,filename))
      for exp in explist:
        exp_ts = exp.imageset.get_image_identifier(0).split('/')[-1].strip()
        idx_filelist.append(exp_ts)
    cbf.append(set(idx_filelist))
  return cbf


def get_indexed_ts_from_cbf(roots):
  ''' Function to get timestamps of images indexed in provided folders. Based on CBF filename '''
  from dxtbx.model.experiment_list import ExperimentListFactory
  cbf = []
  for root in roots:
    idx_filelist = []
    for filename in os.listdir(root):
      if os.path.splitext(filename)[1] != ".cbf": continue
      idx_filelist.append(filename)
    cbf.append(set(idx_filelist))
  return cbf

def plot_venn(params):
  roots = []
  tags = []
  for path in params.input_path:
    roots.append(os.path.abspath(os.path.join(path, 'out')))
    tags.append(path.strip().split('/')[-1])

  results = get_indexed_ts(roots)
  if len(results) == 2:
    try:
      from matplotlib_venn import venn2 as venn_plotter
    except ImportError as e:
      raise Sorry(message)
  elif len(results) == 3:
    try:
      from matplotlib_venn import venn3 as venn_plotter
    except ImportError as e:
      raise Sorry(message)
  else:
    raise Sorry('matplotlib_venn does not currently support plotting anything other than 2 or 3 sets')
  fig_object = plt.figure()
  venn_plotter(results, set_labels = tags)

  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s'%params.pickle_filename, fig_object)
  if params.show_plot:
    plt.show()

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '--h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:], phil_scope=venn_phil_scope)
  plot_venn(params)
