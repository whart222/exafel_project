from __future__ import division
import pandas as pd
import sys, os

from exafel_project.nks.json import to_json as tj
import pandas as pd

class pd_table():
  def __init__(self,merge_root, tag):
    self.t1, self.t2 = self.gen_tables(merge_root,tag)

  @staticmethod
  def nested_dict(d, k0):
    d_list=[]
      
    if isinstance(d,dict):
      for k in d.iterkeys():
        d_list.extend(pd_table.nested_dict(d[k],k))
    else:
      d_list.append({k0:d})
    return d_list

  def gen_tables(self,merge_root,tag):

    xmerge_path = os.path.join(merge_root,tag,"%s_mark0.log"%tag)
    # Argument: Full path to cxi.xmerge log file
    I = tj.xmerge_interpreter(xmerge_path)
    pdb_path = os.path.join(merge_root,tag,"%s_001.pdb"%tag)
    P = tj.pdb_interpreter(pdb_path)
    I.set_cell_and_space_group(P)
    I.summary["Resolution"]="%.1f-%s (%s-%s)"%(P.low_res,I.data["Resolution High"][-1],
                          I.data["Resolution Low"][-1],I.data["Resolution High"][-1])

    xtriage_path = os.path.join(merge_root,tag,"xtriage_%s.out"%tag)
    X = tj.xtriage_interpreter(xtriage_path)
    I.summary["Wilson B factor"]=X.b_factor

    anomalous_path = os.path.join(merge_root,tag,"%s_peakht.log"%tag)
    X = tj.anomalous_interpreter(anomalous_path)
    I.refi["Anomalous peak height"]=X.peak_height
    I.refi.update(P.refi)

    molprobity_path = os.path.join(merge_root,tag,"%s_molprobity.out"%tag)
    M = tj.molprobity_interpreter(molprobity_path)
    I.refi["No. atoms"] = M.n_atoms
    I.refi["B-factors"] = M.b_factors
  
    #from IPython import embed; embed()
    t2 = pd.DataFrame(I.whole['Table 2'])
  
    t1_d = self.nested_dict(I.whole['Table 1'],'Table 1')
    t1_res = {}
    for k in t1_d:
      t1_res.update(k)
    
    t1 = pd.DataFrame(t1_res,index=[0])
    return t1,t2
