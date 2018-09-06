from __future__ import division
from six.moves import range
import sys
from scitbx.simplex import simplex_opt
import scitbx.lbfgs

class lbfgs_helper():
  def __init__(self, x_obs, y_obs, w_obs, initial):
    assert x_obs.size() == y_obs.size()
    self.x_obs = x_obs
    self.y_obs = y_obs
    self.w_obs = w_obs
    self.n = len(initial)
    self.NN = len(initial)//2
    self.x = initial.deep_copy()
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
                     termination_params = scitbx.lbfgs.termination_parameters
                     (traditional_convergence_test_eps=1.e-2, min_iterations=0), core_params =
                     scitbx.lbfgs.core_parameters(gtol=0.1), log=sys.stdout)
    self.a = self.x

  def compute_functional_and_gradients(self):
    self.a = self.x
    f,g = self.target_func_and_grad()
    return f,g

  def target_func_and_grad(self):
    from scitbx.array_family import flex
    import math
    coord_x = self.x[0:self.NN]
    coord_y = self.x[self.NN:2*self.NN]
    outer_coord_x = coord_x.matrix_outer_product(coord_x)
    outer_coord_y = coord_y.matrix_outer_product(coord_y)
    inner = self.y_obs-outer_coord_x-outer_coord_y
    elements = inner*inner
    result = flex.sum(elements)

    grad = -2.*inner.matrix_multiply(coord_x)
    grad = grad.concatenate(-2.*inner.matrix_multiply(coord_y))

    return result,grad
    # ------- Original Implementation :: Super slow ------------------
    # Note that plotting should change since here we assume self.x is [x1,y1,x2,y2, .....]
    #grad = flex.double(self.n)
    #result = 0.0
    #NN = self.x.size()//2
    #for i in range(0,NN-1):
    #  for j in range(i+1,NN):
    #    dx = self.x[i]-self.x[j]
    #    dy = self.x[i+NN]-self.x[j+NN]
    #    d = math.sqrt(dx*dx+dy*dy)
    #    y_diff = 1-self.y_obs[(i*NN)+(j)]-d
    #    grad[i] += -2*y_diff*(self.x[i]-self.x[j])/d
    #    grad[i+NN] += -2*y_diff*(self.x[i+NN]-self.x[j+NN])/d
    #    #grad[i+NN] += -2*y_diff*self.x[j+1]
    #    result += y_diff*y_diff
    #return result,grad
    # -------------- End original implementation ----------------------

class SimplexMinimizer(object):
  def __init__(self, r, x, seed=None, plot=False):
    from dials.array_family import flex
    self.plot=plot
    assert seed is not None, 'seed should not be None'
    self.n = r.focus()[0]*2 # Number of parameters
    self.x = x
    self.r = r
    self.starting_simplex = []
    mt = flex.mersenne_twister(seed)
    random_scale = .5
    self.n_datapoints = r.focus()[0]
    self.data = r
    for i in range(self.n+1):
      self.starting_simplex.append(random_scale*(((mt.random_double(self.n))/2.0)-1.0)+self.x)

    self.optimizer = simplex_opt(dimension=self.n, matrix=self.starting_simplex, evaluator=self, tolerance=1e-3)
    self.x = self.optimizer.get_solution()

  def target(self, vector):
    f2 = 0.0
    NN = self.r.focus()[0]
    for i in range(0,len(vector)-2,2):
      for j in range(i+2, len(vector),2):

        #from IPython import embed; embed(); exit()
        f= self.r[(i//2)*NN+(j//2)] - (vector[i]*vector[j]+vector[i+1]*vector[j+1])
        f2 += f*f
    print ('inside target = %d'%f2)
    return f2


def plot_with_dimensional_embedding(r, show_plot=True, use_lbfgs=True):
  '''
  Plots a high dimensional vector in a 2-d plane using the idea from Diedrichs and Brehms (2014)
  Equation to be minimized is phi = \sum[(rij-xi.xj)^2]
  r = distance metric used. Should be 1-normalized_distance_matrix
  '''
  #
  from scitbx.math import flex
  assert r.focus()[0] == r.focus()[1], 'r matrix has to be square'
  x = []
  flex.set_random_seed(22)

  x = flex.random_double(2*r.focus()[0])
  flex.set_random_seed(22)
  if use_lbfgs:
    print 'Starting LBFGS'
    w_obs = flex.double([1.0]*(r.size()))
    fit = lbfgs_helper(x_obs=r, y_obs=r, w_obs=w_obs,initial=x)
    x = fit.x
    xx = x[0:r.focus()[0]]
    yy = x[r.focus()[0]:2*r.focus()[0]]
  else:
    x = SimplexMinimizer(r,x,seed=22).x
    xx = []
    yy = []
    for i,elem in enumerate(x):
      if i%2 == 0:
        xx.append(elem)
      else:
        yy.append(elem)

  if show_plot:
    import matplotlib.pyplot as plt
    plt.figure(2)
    plt.scatter(xx,yy,c='g',marker='^')
    #plt.xlim([-1,1])
    #plt.ylim([-1,1])
    plt.show()

def run_detail(show_plot, save_plot, use_dummy_data=False):
    file_name = sys.argv[1]
    from xfel.clustering.singleframe import CellOnlyFrame
    from cctbx import crystal
    cells = []
    for line in open(file_name, "r").xreadlines():
      tokens = line.strip().split()
      unit_cell = tuple(float(x) for x in tokens[0:6])
      space_group_symbol = tokens[6]
      crystal_symmetry = crystal.symmetry(unit_cell = unit_cell, space_group_symbol = space_group_symbol)
      cells.append(CellOnlyFrame(crystal_symmetry,path=None))
    MM = [c.mm for c in cells] # get all metrical matrices
    from scitbx.array_family import flex
    MM_double = flex.double()
    for i in range(len(MM)):
      Tup = MM[i]
      for j in range(6):  MM_double.append(Tup[j])

    print("There are %d cells"%(len(MM)))
    if show_plot or save_plot:
      import matplotlib
      if not show_plot:
        # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
        matplotlib.use('Agg') # use a non-interactive backend
      from matplotlib import pyplot as plt
      plt.figure(1)
      plt.plot([c.uc[0] for c in cells],[c.uc[1] for c in cells],"k.", markersize=3.)
      plt.axes().set_aspect("equal")
      if save_plot:
        plt.savefig(plot_name,
                    size_inches=(10,10),
                    dpi=300,
                    bbox_inches='tight')
      if show_plot:
        plt.show()

    print ("Now constructing a Dij matrix.")
    NN = len(MM)
    import omptbx
    omptbx.omp_set_num_threads(64)
    from cctbx.uctbx.determine_unit_cell import NCDist_flatten
    if use_dummy_data:
      '''
      Generate blob data using sklearn. See example here.
      http://scikit-learn.org/stable/auto_examples/cluster/plot_cluster_comparison.html
      '''
      try:
        from sklearn import datasets
      except ImportError:
        print ("Module sklearn not available. Needed to generate dummy data.")
      import numpy as np
      NN = 100
      blobs = datasets.make_blobs(n_samples=NN, random_state=22)
      xx = []
      yy = []
      Dij = np.zeros([NN,NN])
      Dij = flex.double(Dij)
      for x,y in blobs[0]:
        xx.append(x)
        yy.append(y)
      # Get Dij matrix
      for i in range(len(xx)):
        for j in range(len(xx)):
          dij2 = (xx[i]-xx[j])*(xx[i]-xx[j]) + (yy[i]-yy[j])*(yy[i]-yy[j])
          dij = np.sqrt(dij2)
          Dij[i*len(xx)+j] = dij
      if show_plot:
        import matplotlib.pyplot as plt
        #plt.figure()
        plt.scatter(xx,yy)
        plt.show()
    else:
      Dij = NCDist_flatten(MM_double) # loop is flattened
    plot_with_dimensional_embedding(Dij/flex.max(Dij), show_plot=show_plot)

if __name__ == "__main__":
  '''
  Please provide a filename at the command line whose data is to be parsed.
  Example files are given in the exafel_project/ADSE13-25 itself
  '''
  run_detail(show_plot=True, save_plot=False, use_dummy_data=False)
