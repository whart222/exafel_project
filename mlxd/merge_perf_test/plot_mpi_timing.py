import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

start = np.loadtxt(open(str(sys.argv[1])), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')})
end = np.loadtxt(open(str(sys.argv[2])), delimiter=',', dtype={'names': ('A','B','t'), 'formats':('i4','i4','f8')})

ds=[{'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(start['A'],start['B'],start['t']) }]
de=[{'%s:%s'%(a,b): (a,b,t) for a,b,t in zip(end['A'],end['B'],end['t']) }]

#3D Rank A:B vs time diagram
fig = plt.figure()
fig.clf()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlabel('time [t]')
ax.set_ylabel('Rank To Merge')
ax.set_xlabel('Rank Base')
for a in ds[0].keys():
    ax.scatter(ds[0][a][0],ds[0][a][1],ds[0][a][2]-np.min(start['t']),c='r',marker='o') #Plot start
    ax.scatter(de[0][a][0],de[0][a][1],de[0][a][2]-np.min(start['t']),c='b',marker='x') #Plot end
    ax.plot([ds[0][a][0],de[0][a][0] ], [ds[0][a][1],de[0][a][1] ], [ds[0][a][2]-np.min(start['t']), de[0][a][2]-np.min(start['t'])  ], c='k') #Draw line between start and finish
ax.set_zlim3d([0,np.max(end['t']) - np.min(start['t'])])
ax.set_ylim3d([0,60])
ax.set_xlim3d([0,60])
plt.savefig('3d_%s_%s.pdf'%(str(sys.argv[1]),str(sys.argv[2])))
plt.clf()

#2D connections diagram
plt.scatter(start['t'] - np.min(start['t']), start['A'], marker='o', c='b')
plt.scatter(start['t'] - np.min(start['t']), start['B'], marker='o', c='b')
plt.scatter(end['t'] - np.min(start['t']), end['A'], marker='x', c='r')
for ii in xrange(len(end['A'])):
    plt.plot( [ start['t'][ii]-np.min(start['t']), end['t'][ii] - np.min(start['t'])] , [start['B'][ii], end['A'][ii]], c='k')
    if ii <= np.max([end['A'],end['B']]):
        plt.axhline(ii, xmin=0, xmax=1)
plt.xlabel('time [t]')
plt.ylabel('MPI rank')
plt.title('%s_%s'%(str(sys.argv[1]),str(sys.argv[2])))
plt.savefig('2d_%s_%s.pdf'%(str(sys.argv[1]),str(sys.argv[2])))
