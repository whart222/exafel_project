from __future__ import print_function
message = '''
Plots venn diagram of indexing results from multiple folders
This module requires a module called matplotlib_venn. 
Please do pip install matplotlib_venn if you don't have it in your python installation
'''

from matplotlib import pyplot as plt
import os
from libtbx.phil import parse
from libtbx.utils import Sorry
from .indexing_analytics import params_from_phil
try:
  from matplotlib_venn import venn3
except ImportError as e:
  raise Sorry(message)

venn_phil_scope = parse('''
  input_path = None
    .multiple = True
    .type = path
    .help = Path to folders whose indexing results are to be compared and plotted as Venn diagram \
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
        idx_filelist.append(exp.imageset.get_image_identifier(0).split('/')[-1])
    cbf.append(set(idx_filelist))
  return cbf

if __name__ == '__main__'
  if '--help' in sys.argv[1:] or '--h' in sys.argv[1:]:
    print (message)
  params = params_from_phil(sys.argv[1:])
  roots = []
  tags = []
  for path in params.input_path:
    roots.append(os.path.abspath(path))
    tags.append(path.strip().split('/')[-1])

  results = get_indexed_ts(roots)
  fig_object = plt.figure()
  venn3(results, set_labels = tags)

  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s'%params.pickle_filename, fig_object)
  if params.show_plot: 
    plt.show()
