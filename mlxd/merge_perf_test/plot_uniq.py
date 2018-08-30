from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os, sys
from matplotlib import rc

cwd = os.getcwd()
title=sys.argv[1] #Hostname arg1
print "Title=%s"%title
ranks=int(sys.argv[2]) #MPI ranks arg2
print "Ranks=%d"%ranks
dataset=sys.argv[3] #TARs used
print "Dataset=%s"%dataset
dataStrs=sys.argv[4] #CSV string of different data to plot. String to be split, and data loaded as $cwd_$str.dat
print "DataCSV=%s"%dataStrs
t0=sys.argv[5] #t0 offset
print "t0=%s"%t0
tag=sys.argv[6] #t0 offset
print "tag=%s"%tag

#Get the different bar identifiers
with open(tag+'_unique.dat') as fUniqDat:
  uniqStr = fUniqDat.readlines()
uniqStr = [x.strip() for x in uniqStr]
print uniqStr
data = {}
counter = len(uniqStr)
for u in uniqStr:
  try:
    #Load data
    dat = { u : np.loadtxt(cwd + '/' +  tag + '_' + u + '.dat',dtype={'names':('StartEnd','time'), 'formats':('i4','f8',)}, delimiter=',')} #Create dictionary of loaded data; key is the specific set
    dat[u]['time'] -= float(t0) #Start clock at 0
    data.update(dat)
    print "Loaded data %s"%u
  except Exception as e:
    print "Inconsistency in %s data set"%u
    print e
    exit(-1)
  for i in xrange(len(data[u]['StartEnd'])//2):
    plt.plot([data[u]['time'][i],data[u]['time'][i+1]],[1+(i/len(data[u]['StartEnd'])),1+(i/len(data[u]['StartEnd']))])
plt.savefig('test.pdf')

t_start = np.array([]);
t_end = np.array([]);
for kk in uniqStr:
  t_start = np.append(t_start, np.min(data[kk][(data[kk]['StartEnd']==0)]['time']))
  t_end = np.append(t_end, np.max(data[kk][(data[kk]['StartEnd']==1)]['time']))

dat_d = t_end - t_start
dat_dS = np.sort(dat_d)
bins = np.arange(len(dat_d))
fig,ax = plt.subplots(1)
p1 = ax.bar(bins, dat_dS, width=0.35, align='center',color=(0.25,0.5,0.8))

plt.rc('text',usetex=True)
plt.rc('font',family='serif')
titleStr = '\noindent ' + title + '; ' + str(ranks) + ' MPI ranks; runs='  + dataset + ' \\%s'
plt.title(repr(titleStr)%cwd.replace("_","\_").rsplit("/",1)[1], fontsize=8)
plt.ylabel('t [s]')
uniqStrUS = [r'%s'%us.replace('_','\_') for us in uniqStr]

#Split strings to fit onto plot
from IPython import embed; embed()
str_abbrev = [ u[0][0:2] if len(u) is 1 else "".join([ a[:][0:2] for a in u]) for u in [u.split('_') for u in [uniqStr[ii] for ii in np.argsort(dat_d)] ]]

plt.xticks(bins,tuple(str_abbrev), fontsize=6)
#plt.show()
'''
labelstr='r\noindent \begin{eqnarray*} '
for label in uniqStr:
  l_tmp='\textrm{' + label + '} \\ '
  labelstr += l_tmp
labelstr+=' \end{eqnarray*} '
'''

#Construct the labels from the individual measured components
legendStr = r"\noindent \begin{eqnarray*}"
counter=0
for ii in np.argsort(dat_d):
  legendStr += r" \textrm{%s} &=& \textrm{%s} \\ "%(str_abbrev[counter],uniqStrUS[ii])
  counter +=1
legendStr += r" \end{eqnarray*}"
print legendStr

#Plot the legends
props0 = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
ax.text(0.05, 0.95, legendStr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props0)

totalsum = np.sum(dat_d)
bars = ax.patches
labels = [r"%.3f \%%"%(100*a/totalsum) for a in dat_dS]

csum = np.cumsum(dat_dS)
for rect, label in zip(bars, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2, height + csum[-1]*0.02, label, ha='center', va='bottom')

p2 = plt.step(bins, csum, c=(0.45,0.55,0.55), where='mid', linewidth=2)

ax.set_xlim(min(bins)-0.5, max(bins)+0.5)
ax.text(6., max(csum)-max(csum)/len(dat_d), r"\noindent \textbf{Cumulative\\time}", fontsize=12 , color=(0.45,0.55,0.55))

plt.savefig(cwd.replace("_","\_").rsplit("/",1)[1] + '_timing.pdf')
