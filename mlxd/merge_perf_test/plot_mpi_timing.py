import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

StartStr = str(sys.argv[1])
EndStr = str(sys.argv[2])

start = np.loadtxt(open(StartStr), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')})
end = np.loadtxt(open(EndStr), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')})

ds=[{'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(start['A'],start['B'],start['t']) }]
de=[{'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(end['A'],end['B'],end['t']) }]

t0 = np.min(start['t'])

#3D Rank A:B vs time diagram
fig = plt.figure()
fig.clf()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlabel('time [t]')
ax.set_ylabel('Rank To Merge')
ax.set_xlabel('Rank Base')
for a in ds[0].keys():
    ax.scatter( ds[0][a][0], ds[0][a][1], ds[0][a][2]-t0, c='r', marker='o') #Plot start
    ax.scatter( de[0][a][0], de[0][a][1], de[0][a][2]-t0, c='b', marker='x') #Plot end
    ax.plot( [ ds[0][a][0], de[0][a][0] ], [ ds[0][a][1], de[0][a][1] ], [ ds[0][a][2] - t0, de[0][a][2] - t0 ], c='k') #Draw line between start and finish
ax.set_zlim3d([ 0, np.max(end['t']) - t0 ])
ax.set_ylim3d([ np.min([end['A'], end['B']]), np.max([end['A'],end['B']]) ])
ax.set_xlim3d([ np.min([end['A'], end['B']]), np.max([end['A'],end['B']]) ])
plt.savefig('3d_%s_%s.pdf'%(StartStr, EndStr))
plt.clf()

#2D connections diagram
plt.scatter(start['t'] - t0, start['A'], marker='o', c='b')
plt.scatter(start['t'] - t0, start['B'], marker='o', c='b')
plt.scatter(end['t'] - t0, end['A'], marker='x', c='r')
for ii in xrange( len(start['A']) ):
    plt.plot( [ start['t'][ii] - t0, end['t'][ii] - t0] , [start['B'][ii], end['A'][ii]], c='k')
    if ii <= np.max([start['A'],start['B']]):
    	plt.axhline(ii, xmin=0, xmax=1)
plt.xlabel('time [t]')
plt.ylabel('MPI rank')
plt.title('%s_%s'%(StartStr, EndStr))
plt.savefig('2d_%s_%s.pdf'%(StartStr, EndStr))
