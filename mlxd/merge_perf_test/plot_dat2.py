from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib import rc

cwd = os.getcwd()
title=sys.argv[1] #Hostname arg1
ranks=sys.argv[2] #MPI ranks arg2
dataset=sys.argv[3] #TARs used
dataStrs=sys.argv[4] #CSV string of different data to plot. String to be split, and data loaded as $cwd_$str.dat
t0=sys.argv[5] #t0 offset

#Get the different bar identifiers
with open(cwd+'_uniq.dat') as fUniqDat:
  uniqStr = fUniqDat.readlines()
uniqStr = [x.strip() for x in uniqStr]

for u in uniqStr:
  #Load data
  dat = { u : np.loadtxt(cwd + '_' + u + '.dat',dtype={'names':('StartEnd','time'), 'formats':('i4','f8',)}, delimiter=',')} #Create dictionary of loaded data; key is the specific set
  dat[u]['time'] -= float(t0) #Start clock at 0

for i in xrange(len(dat[u]['SE'])//2):
  plt.plot([dat['t'][i],dat['t'][i+1]],[1+(i/len(dat['SE'])),1+(i/len(dat['SE']))])
plt.savefig('test.pdf')

dat_s=np.array([]);
dat_f=np.array([]);

for ii in xrange(len(dat)):
    dat_s = np.append(dat_s,(dat[ii][2] - dat[0][2]))
    dat_f = np.append(dat_f,(dat[ii][3] - dat[0][2]))

dat_d = dat_f - dat_s
bins = np.arange(len(dat_d))
fig,ax = plt.subplots(1)
p1 = ax.bar(bins,dat_d, width=0.35, align='center',color=(0.25,0.5,0.8))

plt.rc('text',usetex=True)
plt.rc('font',family='serif')
titleStr = '\noindent ' + title + '; ' + ranks + ' MPI ranks; runs='  + dataset + ' \\%s'
plt.title(repr(titleStr)%cwd.replace("_","\_").rsplit("/",1)[1], fontsize=8)
plt.ylabel('t [s]')
plt.xticks(bins,tuple(uniqStr))
#plt.show()
'''
labelstr='r\noindent \begin{eqnarray*} '
for label in uniqStr:
  l_tmp='\textrm{' + label + '} \\ '
  labelstr += l_tmp
labelstr+=' \end{eqnarray*} '
'''
#labelstr = r'\noindent \begin{eqnarray*} \textrm{STP} &=& \textrm{SETUP} \\ \textrm{BC} &=& \textrm{BROADCAST} \\ \textrm{SWS} &=& \textrm{SCALER\_WORKER\_SETUP} \\ \textrm{SCW} &=& \textrm{SCALER\_WORKERS} \\ \textrm{GAT} &=& \textrm{GATHER} \\ \textrm{SMA} &=& \textrm{SCALER\_MASTER\_ADD} \\ \textrm{SMF} &=& \textrm{SCALER\_MASTER\_FINALISE} \end{eqnarray*}'
#props0 = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
#ax.text(0.05, 0.95, labelstr, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props0)

totalsum = np.sum(dat_d)
bars = ax.patches
labels = [r"%.3f \%%"%(100*a/totalsum) for a in dat_d]

csum = np.cumsum(dat_d)
for rect, label in zip(bars, labels):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2, height + csum[-1]*0.02, label, ha='center', va='bottom')

p2 = plt.step(bins, csum, c=(0.45,0.55,0.55), where='mid', linewidth=2)

ax.set_xlim(min(bins)-0.5, max(bins)+0.5)
ax.text(4., max(csum)-max(csum)/len(dat_d), r"\noindent \textbf{Cumulative\\time}", fontsize=12 , color=(0.45,0.55,0.55))

plt.savefig(cwd.replace("_","\_").rsplit("/",1)[1] + '_timing.pdf')
