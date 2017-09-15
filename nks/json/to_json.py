from __future__ import division
# based on code from Chuck Yoon, LCLS
import json, os, io, sys

class xmerge_interpreter(object):
  def __init__(self, file_in):
    xmergeFile = file_in
    with open(xmergeFile) as f: self.xmerge = f.readlines()
    self.data = {}
    self.refi = {}
    self.summary = {"CC*": "N/A","No. collected images": "N/A"}
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
      if self.xmerge[i] == '\n':
        continue
      elif 'All' in self.xmerge[i]:
        overall_cc_int = float(self.xmerge[i].split()[2].strip("%"))
        overall_r_split = float(self.xmerge[i].split()[7].strip("%"))
        overall_cc_iso = float(self.xmerge[i].split()[4].strip("%"))
      else:
        resolution_low.append("%.2f"%(float(self.xmerge[i].split()[1])))
        resolution_high.append("%.2f"%(float(self.xmerge[i].split()[3])))
        cc_int.append(float(self.xmerge[i].split()[5].strip("%")))
        r_split.append(float(self.xmerge[i].split()[10].strip("%")))
        cc_iso.append(float(self.xmerge[i].split()[7].strip("%")))

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
      if 'Bin  Resolution Range  Completeness   %    multi> multi> n_meas   <I>    <I/sig(I)>' in line:
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


if __name__=="__main__":
  #libtbx.python json/to_json.py ${MERGE_ROOT}/${TAG}/${TAG}_mark0.log
  # ${MERGE_ROOT}/${TAG}/${TAG}_001.pdb
  """start here
  parse out pdb, combine with xmerge stuff
  get the xmerge unit cell ranges
  parse out the xtriage log
  parse the anomalous peak height
  populate the refinement parameters
  propagate this to all three merge scripts"""
  merge_root = sys.argv[1]
  tag = sys.argv[2]

  xmerge_path = os.path.join(merge_root,tag,"%s_mark0.log"%tag)
  # Argument: Full path to cxi.xmerge log file
  I = xmerge_interpreter(xmerge_path)
  pdb_path = os.path.join(merge_root,tag,"%s_001.pdb"%tag)

  file_out = os.path.join(merge_root,tag,"xmerge.json")
  I.emit(file_out)
