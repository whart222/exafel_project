from __future__ import division
from cctbx.array_family import flex
from libtbx import group_args
#
# List of consensus functions to be implemented
# 1. unit cell
# 2. orientational
# 3. spot clique
#

class clustering_manager(group_args):
  def __init__(self, **kwargs):
    group_args.__init__(self, **kwargs)
    print ('finished Dij, now calculating rho_i and density')
    from xfel.clustering import Rodriguez_Laio_clustering_2014 as RL
    R = RL(distance_matrix = self.Dij, d_c = self.d_c)
    self.rho = rho = R.get_rho()
    ave_rho = flex.mean(rho.as_double())
    NN = self.Dij.focus()[0]
    i_max = flex.max_index(rho)
    delta_i_max = flex.max(flex.double([self.Dij[i_max,j] for j in xrange(NN)]))
    rho_order = flex.sort_permutation(rho, reverse=True)
    rho_order_list = list(rho_order)
    self.delta = delta = R.get_delta(rho_order=rho_order, delta_i_max=delta_i_max)
    cluster_id = flex.int(NN, -1) # -1 means no cluster
    delta_order = flex.sort_permutation(delta, reverse=True)
    N_CLUST = 5               # maximum number of points to be considerd a cluster
    MAX_PERCENTILE_RHO = 0.90 # cluster centers have to be in the top percentile
    n_cluster = 0
    for ic in range(NN):
      # test the density & rho
      item_idx = delta_order[ic] 
      if delta[item_idx] < 0.25*delta[delta_order[0]]: # too low to be a medoid
        continue 
      item_rho_order = rho_order_list.index(item_idx)
      if (1.0*item_rho_order)/NN < MAX_PERCENTILE_RHO:
        cluster_id[item_idx] = n_cluster
        n_cluster +=1
    print ('Found %d clusters'%n_cluster)
    for x in range(NN):
      if cluster_id[x] >= 0:
        print ("XC", x,cluster_id[x], rho[x], delta[x])
    self.cluster_id_maxima = cluster_id.deep_copy()
    R.cluster_assignment(rho_order, cluster_id)
    self.cluster_id_full = cluster_id.deep_copy()
    
    halo = flex.bool(NN,False)
    border = R.get_border( cluster_id = cluster_id )

    for ic in range(n_cluster): #loop thru all border regions; find highest density
      this_border = (cluster_id == ic) & (border==True)
      if this_border.count(True)>0:
        highest_density = flex.max(rho.select(this_border))
        halo_selection = (rho < highest_density) & (this_border==True)
        if halo_selection.count(True)>0:
          cluster_id.set_selected(halo_selection,-1)
        core_selection = (cluster_id == ic) & ~halo_selection
        highest_density = flex.max(rho.select(core_selection))
        too_sparse = core_selection & (rho.as_double() < highest_density/10.) # another heuristic
        if too_sparse.count(True)>0:
          cluster_id.set_selected(too_sparse,-1)
    self.cluster_id_final = cluster_id.deep_copy()





def get_uc_consensus(experiments_list):
  '''
  Uses the Rodriguez Laio 2014 method to do a clustering of the unit cells and then vote for the highest
  consensus unit cell. Input needs to be a list of experiments object. 
  Clustering code taken from github.com/cctbx-xfel/cluster_regression
  Returns an experiment object with crystal unit cell from the cluster with the most points
  '''
  cells = []
  from xfel.clustering.singleframe import CellOnlyFrame
  for experiment in experiments_list:
    assert len(experiment.crystals()) == 1, 'Should have only one crystal model'
    crystal_symmetry = experiment.crystals()[0].get_crystal_symmetry()
    cells.append(CellOnlyFrame(crystal_symmetry))
  MM = [c.mm for c in cells] # metrical matrices
  MM_double = flex.double()
  for i in range(len(MM)):
    Tup = MM[i]
    for j in range(6):
      MM_double.append(Tup[j])
  print('There are %d cells'%len(MM))
  print ('Now constructing a Dij matrix')
  NN = len(MM)
  from cctbx.uctbx.determine_unit_cell import NCDist_flatten
  Dij = NCDist_flatten(MM_double) 
  d_c = 5
  CM = clustering_manager(Dij=Dij, d_c=d_c)
  n_cluster = 1+flex.max(CM.cluster_id_final)
  print (len(cells), ' datapoints have been analyzed')
  print ('%d CLUSTERS'%n_cluster)
  for i in range(n_cluster):
    item = flex.first_index(CM.cluster_id_maxima, i)
    print ('Cluster %d central Unit cell = %d'%(i, item))
    cells[item].crystal_symmetry.show_summary()
  return experiments_list[item].crystals()[0]
   
