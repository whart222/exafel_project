from __future__ import division
from six.moves import range
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.phil import parse
import os
#
# List of consensus functions to be implemented
# 1. unit cell
# 2. orientational
# 3. spot clique
#
clustering_iota_phil_str = '''
clustering {
  Z_delta = 2.0
    .type = float
    .help = cutoff for delta values used in clustering. All mediods have to be above this cutoff in delta values
  d_c = 6.13
    .type = float
    .help = d_c parameter used during clustering by unit cell dimensions
  d_c_ori = 0.13
    .type = float
    .help = d_c parameter used during clustering by orientational matrix A
  max_percentile_rho_uc = 0.95
    .type = float
    .help = rho of a data point needs to be above this value to be considered a mediod during uc clustering
  max_percentile_rho_ori = 0.85
    .type = float
    .help = rho of a data point needs to be above this value to be considered a mediod during orientational clustering
  min_datapts = 5
    .type = int
    .help = Minimum number of datapoints in each cluster to be able to be considered \
            for further indexing
}

'''
clustering_iota_scope = parse(clustering_iota_phil_str)

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
        other = cryst1_ori, fractional_length_tolerance = 50.00,
        unimodular_generator_range=1)
    cryst2_ori_best=cryst2_ori.change_basis(best_similarity_transform)
  except Exception as e:
    #print (str(e))
    cryst2_ori_best = cryst2_ori
  #print 'difference z-score = ', cryst1_ori.difference_Z_score(cryst2_ori_best)
  return cryst1_ori.difference_Z_score(cryst2_ori_best)

def estimate_d_c(Dij):
  ''' Estimate the value of d_c using the assumption that each cluster will be gaussian distributed in it's dij values.
      If we can find out how many of those gaussians are there in the Dij distribution, we can get an estimate of the d_c
      from the standard deviation of the individual gaussians'''
  from scitbx.array_family import flex
  Dij_max=max(Dij.as_1d())
  Dij_min=min(Dij.as_1d())
  # Rounding off to closest multiple of 10
  n_slots=(int(Dij_max)//10+1)*10
  if n_slots == 10:
    return 1.0
  hist_data=flex.histogram(Dij.as_1d(), n_slots=n_slots)
  # Divide the data further into bins and see if there are dead zones with data on either sides. 
  # This will indicate that there are 2+ clusters
  y=hist_data.slots()
  x=hist_data.slot_centers()
  moving_avg_bin=[]
  for i in range(0, n_slots, 10):
    moving_avg_bin.append(flex.mean(flex.double(list(y[i:i+10]))))
  # There has to be one cluster close to 0.0, take that as reference point and find out where the next cluster is
  min_avg = min(moving_avg_bin)*2.0
  d_c=1.0
  for i, avg in enumerate(moving_avg_bin):
    if avg<=min_avg:
      d_c= float(i*10.0)
      break
  return d_c

class clustering_manager(group_args):
  def __init__(self, **kwargs):
    group_args.__init__(self, **kwargs)
    print ('finished Dij, now calculating rho_i and density')
    from xfel.clustering import Rodriguez_Laio_clustering_2014 as RL
    R = RL(distance_matrix = self.Dij, d_c = self.d_c)
    #from clustering.plot_with_dimensional_embedding import plot_with_dimensional_embedding
    #plot_with_dimensional_embedding(1-self.Dij/flex.max(self.Dij), show_plot=True)
    if hasattr(self, 'strategy') is False:
      self.strategy='default'
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
    print ('Z_DELTA = ',self.Z_delta)


    pick_top_solution=False
    rho_stdev = flex.mean_and_variance(rho.as_double()).unweighted_sample_standard_deviation()
    delta_stdev = flex.mean_and_variance(delta).unweighted_sample_standard_deviation()
    if rho_stdev !=0.0 and delta_stdev !=0:
      rho_z=(rho.as_double()-flex.mean(rho.as_double()))/(rho_stdev)
      delta_z=(delta-flex.mean(delta))/(delta_stdev)
    else:
      pick_top_solution=True
      if rho_stdev == 0.0:
        centroids = [flex.first_index(delta,flex.max(delta))]
      elif delta_stdev == 0.0:
        centroids = [flex.first_index(rho,flex.max(rho))]

    significant_delta = []
    significant_rho = []
    # Define strategy to decide cluster center here. Only one should be true
    debug_fix_clustering = True
    if self.strategy=='one_cluster':
        debug_fix_clustering=False
        strategy2=True
    if self.strategy=='strategy_3':
        debug_fix_clustering=False
        strategy3=True
        strategy2=False

    if debug_fix_clustering:
      if not pick_top_solution:
        delta_z_cutoff = min(1.0, max(delta_z))
        rho_z_cutoff = min(1.0, max(rho_z))
        for ic in range(NN):
          # test the density & rho
          if delta_z[ic] >= delta_z_cutoff or delta_z[ic] <= -delta_z_cutoff:
            significant_delta.append(ic)
          if rho_z[ic] >= rho_z_cutoff or rho_z[ic] <= -rho_z_cutoff:
            significant_rho.append(ic)
        if True:
          # Use idea quoted in Rodriguez Laio 2014 paper
          # " Thus, cluster centers are recognized as points for which the value of delta is anomalously large."
          centroid_candidates=list(significant_delta)
          candidate_delta_z=flex.double()
          for ic in centroid_candidates:
            if ic == rho_order[0]:
              delta_z_of_rho_order_0=delta_z[ic]
            candidate_delta_z.append(delta_z[ic])
          i_sorted=flex.sort_permutation(candidate_delta_z, reverse=True)
          # Check that once sorted the top one is not equal to the 2nd or 3rd position
          # If there is a tie, assign centroid to the first one in rho order
          centroids=[]
          # rho_order[0] has to be a centroid
          centroids.append(rho_order[0])
          
          #centroids.append(centroid_candidates[i_sorted[0]])
          for i in range(0, len(i_sorted[:])):
            if centroid_candidates[i_sorted[i]] == rho_order[0]:
              continue
            if delta_z_of_rho_order_0-candidate_delta_z[i_sorted[i]] > 1.0:
              if i>1:
                if -candidate_delta_z[i_sorted[i-1]]+candidate_delta_z[i_sorted[0]] > 1.0:
                  centroids.append(centroid_candidates[i_sorted[i]])
              else:
                centroids.append(centroid_candidates[i_sorted[i]])
            else:
              break
        if False:
          centroid_candidates = list(set(significant_delta).intersection(set(significant_rho)))
          # Now compare the relative orders of the max delta_z and max rho_z to make sure they are within 1 stdev
          centroids = []
          max_delta_z_candidates = -999.9
          max_rho_z_candidates = -999.9
          for ic in centroid_candidates:
            if delta_z[ic] > max_delta_z_candidates:
              max_delta_z_candidates = delta_z[ic]
            if rho_z[ic] > max_rho_z_candidates:
              max_rho_z_candidates = rho_z[ic]
          for ic in centroid_candidates:
            if max_delta_z_candidates - delta_z[ic] < 1.0 and max_rho_z_candidates - rho_z[ic] < 1.0:
              centroids.append(ic)

      #item_idxs = [delta_order[ic] for ic,centroid in enumerate(centroids)]
      item_idxs=centroids
      for item_idx in item_idxs:
        cluster_id[item_idx] = n_cluster
        print ('CLUSTERING_STATS',item_idx,cluster_id[item_idx] )
        n_cluster +=1
        ####
    elif strategy2:
      # Go through list of clusters, see which one has highest joint rank in both rho and delta lists
      # This will only assign one cluster center based on highest product of rho and delta ranks
      product_list_of_ranks=[]
      for ic in range(NN):
        rho_tmp=self.rho[ic]
        delta_tmp=self.delta[ic]
        product_list_of_ranks.append(rho_tmp*delta_tmp)
      import numpy as np
      item_idx=np.argmax(product_list_of_ranks)
      cluster_id[item_idx]=n_cluster # Only cluster assigned
      print ('CLUSTERING_STATS',item_idx,cluster_id[item_idx])
      n_cluster +=1
    elif strategy3:
      # use product of delta and rho and pick out top candidates
      # have to use a significance z_score to filter out the very best
      product_list_of_ranks=flex.double()
      for ic in range(NN):
        rho_tmp=self.rho[ic]
        delta_tmp=self.delta[ic]
        product_list_of_ranks.append(rho_tmp*delta_tmp)
      import numpy as np
      iid_sorted=flex.sort_permutation(product_list_of_ranks, reverse=True)
      cluster_id[iid_sorted[0]]=n_cluster # first point always a cluster
      n_cluster +=1
      print ('CLUSTERING_STATS S3',iid_sorted[0],cluster_id[iid_sorted[0]])
      product_list_of_ranks[iid_sorted[0]]=0.0 # set this to 0.0 so that the mean/stdev does not get biased by one point
      stdev=np.std(product_list_of_ranks)
      mean=np.mean(product_list_of_ranks)
      n_sorted=3
      if stdev == 0.0:
        n_sorted=1
      
      z_critical = 3.0 # 2 sigma significance ?
      # Only go through say 3-4 datapoints 
      # basically there won't be more than 2-3 lattices on an image realistically
      for iid in iid_sorted[1:n_sorted]:
        z_score=(product_list_of_ranks[iid]-mean)/stdev
        if z_score > z_critical:
          cluster_id[iid]=n_cluster
          n_cluster +=1
          print ('CLUSTERING_STATS S3',iid,cluster_id[iid])
        else:
          break # No point going over all points once below threshold z_score

    else:
      for ic in range(NN):
        item_idx = delta_order[ic]
        if ic != 0:
          if delta[item_idx] <= 0.25*delta[delta_order[0]]: # too low to be a medoid
            continue
        item_rho_order = rho_order_list.index(item_idx)
        if (item_rho_order)/NN < MAX_PERCENTILE_RHO:
          cluster_id[item_idx] = n_cluster
          print ('CLUSTERING_STATS',ic,item_idx,item_rho_order,cluster_id[item_idx])
          n_cluster +=1
    ###
#
    print ('Found %d clusters'%n_cluster)
    for x in range(NN):
      if cluster_id[x] >= 0:
        print ("XC", x,cluster_id[x], rho[x], delta[x])
    self.cluster_id_maxima = cluster_id.deep_copy()
    R.cluster_assignment(rho_order, cluster_id, rho)
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

def get_uc_consensus(experiments_list, show_plot=False, save_plot=False, return_only_first_indexed_model=False,finalize_method = 'reindex_with_known_crystal_models', clustering_params = None):
  '''
  Uses the Rodriguez Laio 2014 method to do a hierarchical clustering of the crystal models and
  then vote for the highest consensus crystal mode. Input needs to be a list of experiments object.
  Clustering code taken from github.com/cctbx-xfel/cluster_regression
  Clustering is first done first based on unit cell dimensions. Then for each of the clusters identified,
  a further clustering is done based on orientational matrix A
  '''
  if return_only_first_indexed_model:
    return [experiments_list[0].crystals()[0]], None
  cells = []

  from xfel.clustering.singleframe import CellOnlyFrame
  # Flag for testing Lysozyme data from NKS.Make sure cluster_regression repository is present and configured
  # Program will exit after plots are displayed if this flag is true
  test_nks = False
  if clustering_params is None:
    clustering_params = clustering_iota_scope

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
    clustered_experiments_list = flex.int()
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
  d_c = clustering_params.d_c #five_number_summary(list(Dij))[1]
  d_c = estimate_d_c(Dij)
  #d_c = flex.mean_and_variance(Dij.as_1d()).unweighted_sample_standard_deviation()
  print ('d_c = ',d_c)
  if len(cells) < 5:
    return [experiments_list[0].crystals()[0]], None
  CM = clustering_manager(Dij=Dij, d_c=d_c, max_percentile_rho=clustering_params.max_percentile_rho_uc,Z_delta=clustering_params.Z_delta, strategy='strategy_3')
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
  if do_orientational_clustering:
    print ('IOTA: Starting orientational clustering')
    Dij_ori = {} # dictionary to store Dij for each cluster
    uc_experiments_list = {} # dictionary to store experiments_lists for each cluster
    from collections import Counter
    uc_cluster_count = Counter(list(CM.cluster_id_final))
    # instantiate the Dij_ori flat 1-d array
    # Put all experiments list from same uc cluster together
    if True:
      from scitbx.matrix import sqr
      from cctbx_orientation_ext import crystal_orientation
      crystal_orientation_list = []
      all_A= []
      for i in range(len(experiments_list)):
        crystal_orientation_list.append(crystal_orientation(experiments_list[i].crystals()[0].get_A(), True))
        #exit()
        A_direct = sqr(crystal_orientation_list[i].reciprocal_matrix()).transpose().inverse()
        all_A.append(A_direct[0])
        print ("Direct A matrix 1st element = %12.6f %12.6f %12.6f"%(A_direct[0], A_direct[1], A_direct[2]))
    #  exit()
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
      if uc_cluster_count[cluster] < clustering_params.min_datapts:
        continue
      Dij_ori[cluster] = flex.double([[0.0]*uc_cluster_count[cluster]]*uc_cluster_count[cluster])
    # Now populate the Dij_ori array
      N_samples_in_cluster = len(uc_experiments_list[cluster])
      for i in range(N_samples_in_cluster-1):
        for j in range(i+1, N_samples_in_cluster):
          dij_ori = get_dij_ori(uc_experiments_list[cluster][i].crystals()[0],uc_experiments_list[cluster][j].crystals()[0])
          A_direct_i = sqr(uc_experiments_list[cluster][i].crystals()[0].get_A()).transpose().inverse()
          A_direct_j = sqr(uc_experiments_list[cluster][j].crystals()[0].get_A()).transpose().inverse()
          #print ("Direct A matrix 1st element = %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f"%(dij_ori, A_direct_i[0], A_direct_j[0], A_direct_i[1],A_direct_j[1], A_direct_i[2], A_direct_j[2] ))
          Dij_ori[cluster][N_samples_in_cluster*i+j] = dij_ori
          Dij_ori[cluster][N_samples_in_cluster*j+i] = dij_ori

    # Now do the orientational cluster analysis
    d_c_ori = clustering_params.d_c_ori # 0.13
    from exafel_project.ADSE13_25.clustering.plot_with_dimensional_embedding import plot_with_dimensional_embedding
    #plot_with_dimensional_embedding(1-Dij_ori[1]/flex.max(Dij_ori[1]), show_plot=True)
    A_matrices = []
    for cluster in Dij_ori:
      #if cluster == 2:
      #  CM_ori = clustering_manager(Dij=Dij_ori[cluster], d_c=d_c_ori, max_percentile_rho=0.85, debug=True)
      d_c_ori = estimate_d_c(Dij_ori[cluster])
      #else:
      #d_c_ori=flex.mean_and_variance(Dij_ori[cluster].as_1d()).unweighted_sample_standard_deviation()
      print ('d_c_ori=',d_c_ori)
      CM_ori = clustering_manager(Dij=Dij_ori[cluster], d_c=d_c_ori, max_percentile_rho=clustering_params.max_percentile_rho_ori, Z_delta=clustering_params.Z_delta, strategy='strategy_3')
      n_cluster_ori = 1+flex.max(CM_ori.cluster_id_final)
      #from IPython import embed; embed(); exit()
      for i in range(n_cluster_ori):
        if len([zz for zz in CM_ori.cluster_id_final if zz == i]) < clustering_params.min_datapts:
          continue
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
        print (A_direct)
      if show_plot:
        # Decision graph
        stretch_plot_factor = 1.05 # (1+fraction of limits by which xlim,ylim should be set)
        import matplotlib.pyplot as plt
        plt.plot(CM_ori.rho, CM_ori.delta, "r.", markersize=3.)
        for x in range(len(list(CM_ori.cluster_id_final))):
          if CM_ori.cluster_id_maxima[x] >=0:
            plt.plot([CM_ori.rho[x]], [CM_ori.delta[x]], "ro")
        #exit()
        plt.xlim([-10,stretch_plot_factor*flex.max(CM_ori.rho)])
        plt.ylim([-10,stretch_plot_factor*flex.max(CM_ori.delta)])
        plt.show()
  # FIXME Still to be worked out what exactly should be returned
  #if return_only_first_indexed_model:
  #  return [experiments_list[0].crystals()[0]], clustered_experiments_list
  # Make sure the crystal models are not too close to each other
  # FIXME should be a PHIL
  #from IPython import embed; embed(); exit()
  min_angle = 5.0 # taken from indexer.py
  close_models_list = []
  # Not used really; other fixes have been made to code to figure out outliers
  # Still keeping this in case it it useful later on. 
  if len(dxtbx_crystal_models) > 10000:
    from dials.algorithms.indexing.compare_orientation_matrices import difference_rotation_matrix_axis_angle
    from cctbx_orientation_ext import crystal_orientation
    from dxtbx.model import Crystal
    for i_a in range(0,len(dxtbx_crystal_models)-1):
      for i_b in range(i_a+1,len(dxtbx_crystal_models)):
        cryst_a = dxtbx_crystal_models[i_a]
        cryst_b = dxtbx_crystal_models[i_b]
        cryst_a_ori = crystal_orientation(cryst_a.get_A(), True)
        cryst_b_ori = crystal_orientation(cryst_b.get_A(), True)
        try:
          best_similarity_transform = cryst_b_ori.best_similarity_transformation(
            other = cryst_a_ori, fractional_length_tolerance = 20.00,
            unimodular_generator_range=1)
          cryst_b_ori_best=cryst_b_ori.change_basis(best_similarity_transform)
        except Exception as e:
          cryst_b_ori_best = cryst_b_ori

        # FIXME hardcoded space group for myoglobin LS49
        cryst_b_best=Crystal(cryst_b_ori_best.direct_matrix()[0:3], cryst_b_ori_best.direct_matrix()[3:6], cryst_b_ori_best.direct_matrix()[6:9], 'P 1 21 1')
        R_ab, axis, angle, cb_op_ab = difference_rotation_matrix_axis_angle(cryst_a, cryst_b_best)
        # FIXME
        if abs(angle) < min_angle: # degrees
          close_models_list.append((i_a, i_b))

  # Now prune the dxtbx_crystal_models list
    unique_experiments_list=flex.int(range(len(dxtbx_crystal_models)))
    for close_models in close_models_list:
      i_a,i_b = close_models
      if dxtbx_crystal_models[i_a] is not None and dxtbx_crystal_models[i_b] is not None:
        dxtbx_crystal_models[i_b]=None
        unique_experiments_list[i_b]=i_a
        clustered_experiments_list.set_selected(clustered_experiments_list==i_b, i_a)

    counter=-1
    for ii,model in enumerate(dxtbx_crystal_models):
      if model is not None:
        counter +=1 
        clustered_experiments_list.set_selected(clustered_experiments_list==unique_experiments_list[ii], counter)
    dxtbx_crystal_models=[x for x in dxtbx_crystal_models if x is not None]

  #from IPython import embed; embed(); exit()
  if len(dxtbx_crystal_models) > 0:
    return dxtbx_crystal_models, list(clustered_experiments_list)
  else:
    # If nothing works, atleast return the 1st crystal model that was found
    return [experiments_list[0].crystals()[0]], None
