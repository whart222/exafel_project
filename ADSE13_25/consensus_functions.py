from __future__ import division
from six.moves import range
from cctbx.array_family import flex
from libtbx import group_args
import os
#
# List of consensus functions to be implemented
# 1. unit cell
# 2. orientational
# 3. spot clique
#

def get_dij_ori(cryst1, cryst2, is_reciprocal=True):
  '''
  Takes in 2 dxtbx crystal models, returns the distance between the 2 models in crystal
  space. Currently the distance is defined as the Z-score difference as implemented in
  cctbx/uctbx
  Info regarding some params in best_similarity_transformation after discussion with NKS
  fractional_length_tolerance :: could have value 1 or 200 ??
  unimodular_generator_range = to keep the volume to be 1. If volume doubles on change of basis, make it 2

  '''
  from scitbx.math import flex
  from cctbx_orientation_ext import crystal_orientation
  cryst1_ori = crystal_orientation(cryst1.get_A(), is_reciprocal)
  cryst2_ori = crystal_orientation(cryst2.get_A(), is_reciprocal)
  try:
    best_similarity_transform = cryst2_ori.best_similarity_transformation(
        other = cryst1_ori, fractional_length_tolerance = 10.00,
        unimodular_generator_range=1)
    cryst2_ori_best=cryst2_ori.change_basis(best_similarity_transform)
  except Exception:
    cryst2_ori_best = cryst2_ori
  #print 'difference z-score = ', cryst1_ori.difference_Z_score(cryst2_ori_best)
  return cryst1_ori.difference_Z_score(cryst2_ori_best)


class clustering_manager(group_args):
  def __init__(self, **kwargs):
    group_args.__init__(self, **kwargs)
    print ('finished Dij, now calculating rho_i and density')
    from xfel.clustering import Rodriguez_Laio_clustering_2014 as RL
    R = RL(distance_matrix = self.Dij, d_c = self.d_c)
    #from IPython import embed; embed()
    #from clustering.plot_with_dimensional_embedding import plot_with_dimensional_embedding
    #plot_with_dimensional_embedding(1-self.Dij/flex.max(self.Dij), show_plot=True)
    self.rho = rho = R.get_rho()
    ave_rho = flex.mean(rho.as_double())
    NN = self.Dij.focus()[0]
    i_max = flex.max_index(rho)
    delta_i_max = flex.max(flex.double([self.Dij[i_max,j] for j in range(NN)]))
    rho_order = flex.sort_permutation(rho, reverse=True)
    rho_order_list = list(rho_order)
    self.delta = delta = R.get_delta(rho_order=rho_order, delta_i_max=delta_i_max)
    cluster_id = flex.int(NN, -1) # -1 means no cluster
    delta_order = flex.sort_permutation(delta, reverse=True)
    MAX_PERCENTILE_RHO = self.max_percentile_rho # cluster centers have to be in the top percentile
    n_cluster = 0
#
#
    for ic in range(NN):
      # test the density & rho
      item_idx = delta_order[ic]
      delta_stats = flex.mean_and_variance(delta)
      Z_delta = 2.0
      if ic != 0:
        if (delta[item_idx]-delta_stats.mean())/delta_stats.unweighted_sample_standard_deviation() < Z_delta:
        #if delta[item_idx] <= 0.25*delta[delta_order[0]]: # too low to be a medoid
          continue
      try:
        if self.debug:
          from IPython import embed; embed()
      except AttributeError:
        print ('debug attribute not present')
      item_rho_order = rho_order_list.index(item_idx)
      if (item_rho_order)/NN < MAX_PERCENTILE_RHO:
        cluster_id[item_idx] = n_cluster
        print ('CLUSTERING_STATS',ic,item_idx,item_rho_order,cluster_id[item_idx])
        n_cluster +=1
#
#
    print ('Found %d clusters'%n_cluster)
    for x in range(NN):
      if cluster_id[x] >= 0:
        print ("XC", x,cluster_id[x], rho[x], delta[x])
    self.cluster_id_maxima = cluster_id.deep_copy()
    #import pdb; pdb.set_trace()
    R.cluster_assignment(rho_order, cluster_id)
    self.cluster_id_full = cluster_id.deep_copy()

    #halo = flex.bool(NN,False)
    #border = R.get_border( cluster_id = cluster_id )

    #for ic in range(n_cluster): #loop thru all border regions; find highest density
    #  this_border = (cluster_id == ic) & (border==True)
    #  if this_border.count(True)>0:
    #    highest_density = flex.max(rho.select(this_border))
    #    halo_selection = (rho < highest_density) & (this_border==True)
    #    if halo_selection.count(True)>0:
    #      cluster_id.set_selected(halo_selection,-1)
    #    core_selection = (cluster_id == ic) & ~halo_selection
    #    highest_density = flex.max(rho.select(core_selection))
    #    too_sparse = core_selection & (rho.as_double() < highest_density/10.) # another heuristic
    #    if too_sparse.count(True)>0:
    #      cluster_id.set_selected(too_sparse,-1)
    self.cluster_id_final = cluster_id.deep_copy()

def get_uc_consensus(experiments_list, show_plot=False, save_plot=False, return_only_first_indexed_model=False,finalize_method = 'reindex_with_known_crystal_models'):
  '''
  Uses the Rodriguez Laio 2014 method to do a clustering of the unit cells and then vote for the highest
  consensus unit cell. Input needs to be a list of experiments object.
  Clustering code taken from github.com/cctbx-xfel/cluster_regression
  '''
  cells = []

  from xfel.clustering.singleframe import CellOnlyFrame
  # Flag for testing Lysozyme data from NKS.Make sure cluster_regression repository is present and configured
  # Program will exit after plots are displayed if this flag is true
  test_nks = False
  if test_nks:
    from cctbx import crystal
    import libtbx.load_env
    cluster_regression = libtbx.env.find_in_repositories(
        relative_path="cluster_regression",
        test=os.path.isdir)
    file_name = os.path.join(cluster_regression, 'examples', 'lysozyme1341.txt')
    for line in open(file_name, "r").xreadlines():
      tokens = line.strip().split()
      unit_cell = tuple(float(x) for x in tokens[0:6])
      space_group_symbol = tokens[6]
      crystal_symmetry = crystal.symmetry(unit_cell = unit_cell, space_group_symbol = space_group_symbol)
      cells.append(CellOnlyFrame(crystal_symmetry))
  else:
    clustered_experiments_list = []
    for experiment in experiments_list:
      if len(experiment.crystals()) >1: print ('IOTA:Should have only one crystal model')
      crystal_symmetry = experiment.crystals()[0].get_crystal_symmetry()
      cells.append(CellOnlyFrame(crystal_symmetry))
      # Maintain a list which is meaningless right now that will finally contain the
      # final clustering results
      clustered_experiments_list.append(-1)
  MM = [c.mm for c in cells] # metrical matrices
  MM_double = flex.double()
  for i in range(len(MM)):
    Tup = MM[i]
    for j in range(6):
      MM_double.append(Tup[j])
  print('There are %d cells'%len(MM))
  coord_x = flex.double([c.uc[0] for c in cells])
  coord_y = flex.double([c.uc[1] for c in cells])
  if show_plot or save_plot:
    import matplotlib
    if not show_plot:
      matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.plot([c.uc[0] for c in cells],[c.uc[1] for c in cells],"k.", markersize=3.)
    plt.axes().set_aspect("equal")
  if save_plot:
    plot_name = 'uc_cluster.png'
    plt.savefig(plot_name,
                size_inches=(10,10),
                dpi=300,
                bbox_inches='tight')
  if show_plot:
    plt.show()
  print ('Now constructing a Dij matrix: Starting Unit Cell clustering')
  NN = len(MM)
  from cctbx.uctbx.determine_unit_cell import NCDist_flatten
  Dij = NCDist_flatten(MM_double)
  from scitbx.math import five_number_summary
  d_c = 6.13 #five_number_summary(list(Dij))[1]
  CM = clustering_manager(Dij=Dij, d_c=d_c, max_percentile_rho=0.95)
  n_cluster = 1+flex.max(CM.cluster_id_final)
  print (len(cells), ' datapoints have been analyzed')
  print ('%d CLUSTERS'%n_cluster)
  for i in range(n_cluster):
    item = flex.first_index(CM.cluster_id_maxima, i)
    print ('Cluster %d central Unit cell = %d'%(i, item))
    cells[item].crystal_symmetry.show_summary()

  # More plots for debugging
  appcolors = ['b', 'r', '#ff7f0e', '#2ca02c',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
  if show_plot:
    # Decision graph
    import matplotlib.pyplot as plt
    #from IPython import embed; embed()
    plt.plot(CM.rho, CM.delta, "r.", markersize=3.)
    for x in range(NN):
      if CM.cluster_id_maxima[x] >=0:
        plt.plot([CM.rho[x]], [CM.delta[x]], "ro")
    plt.show()

  if show_plot:
    import matplotlib.pyplot as plt
    colors = [appcolors[i%10] for i in CM.cluster_id_full]
    plt.scatter(coord_x, coord_y, marker='o', color=colors, linewidth=0.4, edgecolor='k')
    for i in range(n_cluster):
      item = flex.first_index(CM.cluster_id_maxima, i)
      plt.plot([cells[item].uc[0]], cells[item].uc[1], 'y.')
      plt.axes().set_aspect("equal")
      plt.show()
  if test_nks:
    exit()

  # Now look at each unit cell cluster for orientational clustering
  # idea is to cluster the orientational component in each of the unit cell clusters
  #
  do_orientational_clustering = not return_only_first_indexed_model # temporary.
  dxtbx_crystal_models = []
  #from IPython import embed; embed()
  if do_orientational_clustering:
    print ('IOTA: Starting orientational clustering')
    Dij_ori = {} # dictionary to store Dij for each cluster
    uc_experiments_list = {} # dictionary to store experiments_lists for each cluster
    from collections import Counter
    uc_cluster_count = Counter(list(CM.cluster_id_final))
    # instantiate the Dij_ori flat 1-d array
    #from IPython import embed; embed();
    # Put all experiments list from same uc cluster together
    if True:
      from scitbx.matrix import sqr
      from cctbx_orientation_ext import crystal_orientation
      #crystal_orientation_list = []
      #for i in range(len(experiments_list)):
      #  crystal_orientation_list.append(crystal_orientation(experiments_list[i].crystals()[0].get_A(), True))
        #from IPython import embed; embed();
        #exit()
        #A_direct = sqr(crystal_orientation_list[i].reciprocal_matrix()).transpose().inverse()
        #print ("Direct A matrix 1st element = %12.6f"%A_direct[0])
    CM_mapping = {}
    for i in range(len(experiments_list)):
      if CM.cluster_id_full[i] not in uc_experiments_list:
        uc_experiments_list[CM.cluster_id_full[i]] = []
        CM_mapping[CM.cluster_id_full[i]] = []
      uc_experiments_list[CM.cluster_id_full[i]].append(experiments_list[i])
      # Maintain mapping between original experiments_list and uc_exeriments_list
      # Mapping: key> index_in_experiments_list | value> cluster_id, index_in_uc_cluster
      CM_mapping[CM.cluster_id_full[i]].append((i,len(uc_experiments_list[CM.cluster_id_full[i]])-1))
    for cluster in uc_cluster_count:
      # Make sure there are atleast a minimum number of samples in the cluster
      if uc_cluster_count[cluster] < 5:
        continue
      Dij_ori[cluster] = flex.double([[0.0]*uc_cluster_count[cluster]]*uc_cluster_count[cluster])
    # Now populate the Dij_ori array
      N_samples_in_cluster = len(uc_experiments_list[cluster])
      #from IPython import embed; embed();
      for i in range(N_samples_in_cluster-1):
        for j in range(i+1, N_samples_in_cluster):
          dij_ori = get_dij_ori(uc_experiments_list[cluster][i].crystals()[0],uc_experiments_list[cluster][j].crystals()[0])
          Dij_ori[cluster][N_samples_in_cluster*i+j] = dij_ori
          Dij_ori[cluster][N_samples_in_cluster*j+i] = dij_ori

    # Now do the orientational cluster analysis
    d_c_ori = 0.13
    from exafel_project.ADSE13_25.clustering.plot_with_dimensional_embedding import plot_with_dimensional_embedding
    #plot_with_dimensional_embedding(1-Dij_ori[1]/flex.max(Dij_ori[1]), show_plot=True)
    A_matrices = []
    for cluster in Dij_ori:
      #if cluster == 2:
      #  CM_ori = clustering_manager(Dij=Dij_ori[cluster], d_c=d_c_ori, max_percentile_rho=0.85, debug=True)
      #else:
      CM_ori = clustering_manager(Dij=Dij_ori[cluster], d_c=d_c_ori, max_percentile_rho=0.85)
      n_cluster_ori = 1+flex.max(CM_ori.cluster_id_final)
      for i in range(n_cluster_ori):
        item = flex.first_index(CM_ori.cluster_id_maxima, i)
        dxtbx_crystal_model = uc_experiments_list[cluster][item].crystals()[0]
        dxtbx_crystal_models.append(dxtbx_crystal_model)
        # Map the orientational clusters to the original experiments_list indices
        # This should be the final list of clusters!
        for j,ori_cluster_id in enumerate(CM_ori.cluster_id_final):
          if ori_cluster_id == i:
            xx,yy = CM_mapping[cluster][j]
            clustered_experiments_list[xx] = len(dxtbx_crystal_models)-1
        from scitbx.matrix import sqr
        from cctbx_orientation_ext import crystal_orientation
        crystal_orientation = crystal_orientation(dxtbx_crystal_model.get_A(), True)
        A_direct = sqr(crystal_orientation.reciprocal_matrix()).transpose().inverse()
        A_matrices.append(A_direct)
        print ("IOTA: Direct A matrix 1st element of orientational cluster %d  = %12.6f"%(i,A_direct[0]))
      if show_plot:
        # Decision graph
        stretch_plot_factor = 1.05 # (1+fraction of limits by which xlim,ylim should be set)
        import matplotlib.pyplot as plt
        plt.plot(CM_ori.rho, CM_ori.delta, "r.", markersize=3.)
        for x in range(len(list(CM_ori.cluster_id_final))):
          if CM_ori.cluster_id_maxima[x] >=0:
            plt.plot([CM_ori.rho[x]], [CM_ori.delta[x]], "ro")
        #from IPython import embed; embed();
        #exit()
        plt.xlim([-10,stretch_plot_factor*flex.max(CM_ori.rho)])
        plt.ylim([-10,stretch_plot_factor*flex.max(CM_ori.delta)])
        plt.show()
  #from IPython import embed; embed()
  # FIXME Still to be worked out what exactly should be returned
  if return_only_first_indexed_model:
    return [experiments_list[0].crystals()[0]], clustered_experiments_list
  if len(dxtbx_crystal_models) > 0:
    return dxtbx_crystal_models, clustered_experiments_list
  else:
    # If nothing works, atleast return the 1st crystal model that was found
    return [experiments_list[0].crystals()[0]], clustered_experiments_list
