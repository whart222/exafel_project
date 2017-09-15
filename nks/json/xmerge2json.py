from __future__ import division
#code is courtesy Chuck Yoon, LCLS
import json, io, sys

# Argument: Full path to cxi.xmerge log file
xmergeFile = sys.argv[1]

with open(xmergeFile) as f: xmerge = f.readlines()

startNum = 0
endNum = 0
for lineNum, line in enumerate(xmerge):
    if 'Table of Scaling Results' in line:
        startNum = lineNum + 6 # Numbers start 6 lines below
    if 'All' in line:
        endNum = lineNum

resolution_low = []
resolution_high = []
cc_int = []
n_int = []
r_split = []
for i in range(startNum, endNum+1):
    if xmerge[i] == '\n':
        continue
    elif 'All' in xmerge[i]:
        resolution_low.append('')
        resolution_high.append('')
        cc_int.append(xmerge[i].split()[2])
        n_int.append(int(xmerge[i].split()[3]))
        r_split.append(xmerge[i].split()[7])
    else:
        resolution_low.append(float(xmerge[i].split()[1]))
        resolution_high.append(float(xmerge[i].split()[3]))
        cc_int.append(xmerge[i].split()[5])
        n_int.append(int(xmerge[i].split()[6]))
        r_split.append(xmerge[i].split()[10])

try:
    to_unicode = unicode
except NameError:
    to_unicode = str

data = {"resolution range low": resolution_low,
       "resolution range high": resolution_high,
       "CCint": cc_int,
       "Nint": n_int,
       "Rsplit": r_split
       }

# Write JSON file
with io.open('xmerge.json', 'w', encoding='utf8') as outfile:
    str_ = json.dumps(data,
                      indent=4, sort_keys=True,
                      separators=(',', ': '), ensure_ascii=False)
    outfile.write(to_unicode(str_))
