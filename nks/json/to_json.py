from __future__ import division
# based on code from Chuck Yoon, LCLS
import json, os, io, sys

class xmerge_interpreter(object):
  def __init__(self, file_in):
    xmergeFile = file_in
    with open(xmergeFile) as f: self.xmerge = f.readlines()
    self.data = {}
    self.refi = {}
    self.summary = {"CC*": "N/A","CCano": "N/A", "No. collected images": "N/A"}
    self.whole = {"Table 1":{"Data collection":self.summary,"Refinement":self.refi},
                  "Table 2":self.data}
    self.get_cchalf_table()
    self.get_lattices()
    self.get_nmeas()

  def get_cchalf_table(self):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.xmerge):
      if 'Table of Scaling Results' in line:
          startNum = lineNum + 6 # Numbers start 6 lines below
      if 'All' in line:
          endNum = lineNum

    resolution_low = []
    resolution_high = []
    cc_int = []
    r_split = []
    cc_iso = []
    for i in range(startNum, endNum+1):
      st = self.xmerge[i]
      if st == '\n':
        continue
      elif 'All' in self.xmerge[i]:
        overall_cc_int = float(st.split()[2].strip("%"))
        overall_r_split = float(st.split()[7].strip("%"))
        overall_cc_iso = float(st.split()[4].strip("%"))
      else:
        resolution_low.append("%.2f"%(float(st.split()[1])))
        resolution_high.append("%.2f"%(float(st.split()[3])))
	#Parsing after the ] due to potential whitespace issues with low data numbers  
        cc_int.append(float(st.split("]")[1].split()[0].strip("%")))
        r_split.append(float(st.split("]")[1].split()[5].strip("%")))
        cc_iso.append(float(st.split("]")[1].split()[2].strip("%")))

    # Lowest resolution possibly set to infinity
    for ix in xrange(len(resolution_low)):
        if float(resolution_low[ix])==-1:  resolution_low[ix]="inf"

    self.data.update( dict([("Resolution Low",resolution_low),
       ("Resolution High",resolution_high),
       ("CC1/2", cc_int),("CCiso", cc_iso),
       ("Rsplit", r_split)
       ])
       )
    self.summary["CCiso"] = "%.1f (%.1f)"%(overall_cc_iso,cc_iso[-1])
    self.summary["CC1/2"] = "%.1f (%.1f)"%(overall_cc_int,cc_int[-1])
    self.summary["Rsplit"] = "%.1f (%.1f)"%(overall_r_split,r_split[-1])

  def get_lattices(self):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.xmerge):
      if 'Bin  Resolution Range # images %accept' in line:
          startNum = lineNum + 2
      if startNum > 0 and 'All' in line:
          endNum = lineNum
          break

    N_images = []

    for i in range(startNum, endNum+1):
      if self.xmerge[i] == '\n':
        continue
      elif 'All' in self.xmerge[i]:
        all_images = (int(self.xmerge[i].split()[1]))
      else:
        N_images.append((int(self.xmerge[i].split()[4])))

    self.data["No. Lattices"] = N_images
    self.summary["No. lattices merged"] = all_images
    self.summary["No. images used"] = all_images

  def get_nmeas(self):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.xmerge):
      if "".join('Bin  Resolution Range  Completeness   %    multi> multi> n_meas   <I>    <I/sig(I)>'.split()) \
      in "".join(line.strip().split()):
          startNum = lineNum + 2
    for iline in xrange(startNum,len(self.xmerge)):
      if 'All' in self.xmerge[iline]:
          endNum = iline
          break

    n_meas = []
    n_unique = []
    n_completeness = []
    f_multiplicity = []
    f_ioversigma = []
    for i in range(startNum, endNum+1):
      if self.xmerge[i] == '\n':
        continue
      elif 'All' in self.xmerge[i]:
        all_n_meas = (int(self.xmerge[i].split()[5]))
        all_completeness = float(self.xmerge[i].split()[2])
        all_multiplicity = float(self.xmerge[i].split()[3])
        all_ioversigma = float(self.xmerge[i].split()[7])
      else:
        n_meas.append((int(self.xmerge[i].split()[8])))
        complete = self.xmerge[i].split()[4][1:-1].split("/")[0]
        n_unique.append(complete)
        n_completeness.append(float(self.xmerge[i].split()[5]))
        f_multiplicity.append(float(self.xmerge[i].split()[6]))
        f_ioversigma.append(float(self.xmerge[i].split()[10]))

    self.data["No. Measurements"] = n_meas
    self.data["No. Unique reflections"] = n_unique
    self.data["Completeness"] = n_completeness
    self.summary["Completeness"] = "%.1f (%.1f)"%(all_completeness,n_completeness[-1])
    self.data["<Multiplicity>"] = f_multiplicity
    self.summary["Multiplicity (Stills)"] = "%.1f (%.1f)"%(all_multiplicity,f_multiplicity[-1])
    self.data["<I/sigI>"] = f_ioversigma
    self.summary["I/sigI"] = "%.1f (%.1f)"%(all_ioversigma,f_ioversigma[-1])
    tot_unique = sum([int(d) for d in n_unique])
    self.summary["No. total reflections"] = "%d (%d)"%(tot_unique,int(n_unique[-1]))

  def set_cell_and_space_group(self,pdb_int):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.xmerge):
      if 'a edge' in line:
          startNum = lineNum + 2
    a_description = "".join(self.xmerge[startNum].split()[1:4])
    for lineNum, line in enumerate(self.xmerge):
      if 'b edge' in line:
          startNum = lineNum + 2
    b_description = "".join(self.xmerge[startNum].split()[1:4])
    for lineNum, line in enumerate(self.xmerge):
      if 'c edge' in line:
          startNum = lineNum + 2
    c_description = "".join(self.xmerge[startNum].split()[1:4])

    self.summary.update(
      {"Space group":pdb_int.space_group,
       "Cell dimensions":{"a":a_description, "b":b_description, "c":c_description,
       "alpha":pdb_int.alpha, "beta":pdb_int.beta, "gamma":pdb_int.gamma}}
    )

  def emit(self,out):
    try:
      self.to_unicode = unicode
    except NameError:
      self.to_unicode = str

    # Write JSON file
    with io.open(out, 'w', encoding='utf8') as outfile:

      str_ = json.dumps(self.whole,
                      indent=4, sort_keys=True,
                      separators=(',', ': '), ensure_ascii=False)
      outfile.write(self.to_unicode(str_))

class pdb_interpreter(object):
  def __init__(self, file_in):
    PdbFile = file_in
    with open(PdbFile) as f: self.pdb = f.readlines()
    self.refi = {}
    self.summary = {}
    self.get_cryst1()
    self.get_other_stuff()

  def get_cryst1(self):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.pdb):
      if 'CRYST1' in line:
          startNum = lineNum
          tokens = line.split()
          print tokens
          self.alpha,self.beta,self.gamma = [float(a) for a in tokens[4:7]]
          self.space_group = " ".join(tokens[7:])
      if 'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS)' in line:
          startNum = lineNum
          tokens = line.strip().split()
          self.low_res = float(tokens[-1])
  def get_other_stuff(self):
    for lineNum, line in enumerate(self.pdb):
      if 'RESOLUTION RANGE HIGH (ANGSTROMS)' in line:
          tokens = line.split()
          self.refi["Resolution"] = float(tokens[-1])
      if 'R VALUE            (WORKING SET)' in line:
          tokens = line.split()
          r_work = float(tokens[-1])*100
      if 'FREE R VALUE                    ' in line:
          tokens = line.split()
          r_free = float(tokens[-1])*100
      if 'REMARK Final:' in line:
          tokens = line.split()
          self.refi["R.m.s deviations"] = {"Bond lengths":float(tokens[10]), "Bond angles": float(tokens[13])}
      if 'REMARK   3    ALL-ATOM CLASHSCORE' in line:
          tokens = line.split()
          self.refi["Clashscore"] = float(tokens[5])
      if 'REMARK   3    RAMACHANDRAN PLOT' in line:
          self.refi["Ramachandran statistics"] = {"Favored":float(self.pdb[lineNum+3].split()[4]),
                                                 "Outliers":float(self.pdb[lineNum+1].split()[4])}
      if 'REMARK         end:' in line:
          tokens = line.split()
          self.refi["No. atoms"] = {"Water":int(tokens[9])}
          self.refi["B-factors"] = {"Water":float(tokens[8])}
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.pdb):
      if 'FIT TO DATA USED IN REFINEMENT (IN BINS)' in line:
          startNum = lineNum + 2
      if startNum > 0 and len(line.split())<4:
          endNum = lineNum - 1
          break
    tokens = self.pdb[endNum].split()
    print tokens
    last_r_work = float(tokens[9])*100
    last_r_free = float(tokens[10])*100
    self.refi["Rwork / Rfree"] = "%.1f/%.1f (%.1f/%.1f)"%(r_work, r_free, last_r_work, last_r_free)

    """
            "No. atoms":{"Protein":5094, "Ligand/ion":5, "Water":452},
            "B-factors":{"Protein":14.7, "Ligand/ion":14.6, "Water":24.1},
     """
class molprobity_interpreter(object):
  def __init__(self, file_in):
    File = file_in
    self.n_atoms = {}
    self.b_factors = {}
    with open(File) as f: self.lines = f.readlines()
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.lines):
      if '  Macromolecules:' in line:
        self.n_atoms["Protein"] = int(self.lines[lineNum+1].split()[4])
        self.b_factors["Protein"] = float(self.lines[lineNum+2].split()[3])
      if '  Ligands:' in line:
        self.n_atoms["Ligand/ion"] = int(self.lines[lineNum+1].split()[4])
        self.b_factors["Ligand/ion"] = float(self.lines[lineNum+2].split()[3])
      if '  Waters:' in line:
        self.n_atoms["Water"] = int(self.lines[lineNum+1].split()[4])
        self.b_factors["Water"] = float(self.lines[lineNum+2].split()[3])

class xtriage_interpreter(object):
  def __init__(self, file_in):
    File = file_in
    with open(File) as f: self.lines = f.readlines()
    self.get_wilson()

  def get_wilson(self):
    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.lines):
      if ' ML estimate of overall B value' in line:
          startNum = lineNum + 1
          break
    self.b_factor = float(self.lines[startNum].split()[0])

class anomalous_interpreter(object):
  def __init__(self, file_in):
    File = file_in
    with open(File) as f: self.lines = f.readlines()

    startNum = 0
    endNum = 0
    for lineNum, line in enumerate(self.lines):
      if 'Maximum grid point value' in line:
          startNum = lineNum + 2
          break
    self.peak_height = float(self.lines[startNum].split()[-2])

if __name__=="__main__":
  #libtbx.python json/to_json.py ${MERGE_ROOT} ${TAG}
  merge_root = sys.argv[1]
  tag = sys.argv[2]

  xmerge_path = os.path.join(merge_root,tag,"%s_mark0.log"%tag)
  # Argument: Full path to cxi.xmerge log file
  I = xmerge_interpreter(xmerge_path)
  pdb_path = os.path.join(merge_root,tag,"%s_001.pdb"%tag)
  P = pdb_interpreter(pdb_path)
  I.set_cell_and_space_group(P)
  I.summary["Resolution"]="%.1f-%s (%s-%s)"%(P.low_res,I.data["Resolution High"][-1],
                          I.data["Resolution Low"][-1],I.data["Resolution High"][-1])

  xtriage_path = os.path.join(merge_root,tag,"xtriage_%s.out"%tag)
  X = xtriage_interpreter(xtriage_path)
  I.summary["Wilson B factor"]=X.b_factor

  anomalous_path = os.path.join(merge_root,tag,"%s_peakht.log"%tag)
  X = anomalous_interpreter(anomalous_path)
  I.refi["Anomalous peak height"]=X.peak_height
  I.refi.update(P.refi)

  molprobity_path = os.path.join(merge_root,tag,"%s_molprobity.out"%tag)
  M = molprobity_interpreter(molprobity_path)
  I.refi["No. atoms"] = M.n_atoms
  I.refi["B-factors"] = M.b_factors

  file_out = os.path.join(merge_root,tag,"%s.json"%tag)
  I.emit(file_out)
