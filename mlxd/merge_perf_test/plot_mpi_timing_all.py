from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.markers as mkr
from matplotlib import colors as mc

StartStr=[]
EndStr=[]
start=[]
end=[]
ds=[]
de=[]

for ii in xrange(len(sys.argv)//2):
  StartStr.append(str(sys.argv[2*ii +1]))
  EndStr.append(str(sys.argv[2*ii +2]))
  start.append( (np.loadtxt( open( StartStr[ii] ), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')}), (StartStr[ii], EndStr[ii]) ))
  end.append( (np.loadtxt( open( EndStr[ii] ), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')}), (StartStr[ii], EndStr[ii]) ))

  ds.append( ({'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(start[ii][0]['A'],start[ii][0]['B'],start[ii][0]['t']) }, (StartStr[ii], EndStr[ii]) ) )
  de.append( ({'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(end[ii][0]['A'],end[ii][0]['B'],end[ii][0]['t']) }, (StartStr[ii], EndStr[ii]) ) )

min_t0=[]
max_t0=[]
min_s=[]
min_e=[]
max_s=[]
max_e=[]
for ii in xrange(len(start)):
  max_t0.append( np.max(start[ii][0]['t']) )
  min_t0.append( np.min(start[ii][0]['t']) )
  min_s.append( np.min( [ start[ii][0]['A'], start[ii][0]['B'] ] ) )
  min_e.append( np.min( [ end[ii][0]['A'], end[ii][0]['B'] ] ) )
  max_s.append( np.max( [ start[ii][0]['A'], start[ii][0]['B'] ] ) )
  max_e.append( np.max( [ end[ii][0]['A'], end[ii][0]['B'] ] ) )

t0 = np.min(min_t0)

#3D Rank A:B vs time diagram
fig = plt.figure()
fig.clf()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlabel('time [s]')
ax.set_ylabel('Rank To Merge')
ax.set_xlabel('Rank Base')
for b in xrange(len(ds)):
  for a in ds[b][0].keys():
    ax.scatter( ds[b][0][a][0], ds[b][0][a][1], ds[b][0][a][2]-t0, c=mc.XKCD_COLORS[mc.XKCD_COLORS.keys()[b]], marker='.', markersize=1, alpha=0.33 )#mkr.MarkerStyle(b), alpha=0.25) #Plot start
    ax.scatter( de[b][0][a][0], de[b][0][a][1], de[b][0][a][2]-t0, c=mc.XKCD_COLORS[mc.XKCD_COLORS.keys()[b*2]], marker='.', markersize=1, alpha=0.33 )#mkr.MarkerStyle(b), alpha=0.25) #Plot end
    ax.plot( [ ds[b][0][a][0], de[b][0][a][0] ], [ ds[b][0][a][1], de[b][0][a][1] ], [ ds[b][0][a][2] - t0, de[b][0][a][2] - t0 ], c='k') #Draw line between start and finish

ax.set_zlim3d([ 0, np.max(max_t0) - t0 ])
ax.set_ylim3d([ np.min([ min_s, min_e ]), np.max([ max_s, max_e ]) ])
ax.set_xlim3d([ np.min([ min_s, min_e ]), np.max([ max_e, max_e ]) ])
plt.savefig('3d_all.pdf')
#plt.savefig('3d_all_%s_%s.pdf'%(StartStr, EndStr))
plt.clf()
exit()
#2D connections diagram
#Draw lines to mark the MPI ranks
for ii in xrange(np.max([start['A'],start['B']])):
  plt.axhline(ii, xmin=0, xmax=1, linewidth=0.5)

#Draw lines between the start and end for reducing 2 data sets
for a in ds[0].keys():

  plt.plot( [ ds[0][a][2] - t0, de[0][a][2] - t0] , [ds[0][a][1], de[0][a][0]], linestyle='-', linewidth=0.5, c='k', alpha=0.8)
  plt.scatter( start['t'] - t0, start['B'], marker='x', c='r', alpha=0.8)
  plt.scatter( end['t'] - t0, end['A'], marker='o', c='b', alpha=0.8)

plt.xlabel('time [s]')
plt.ylabel('MPI rank')
plt.title('%s_%s'%(StartStr, EndStr))
plt.xlim([ 0, np.max(end['t']) - t0 ])
plt.ylim([ np.min([end['A'], end['B']]), np.max([end['A'],end['B']]) ])
plt.savefig('2d_%s_%s.pdf'%(StartStr, EndStr))
