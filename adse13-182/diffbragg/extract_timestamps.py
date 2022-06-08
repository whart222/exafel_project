import os
import sys
#import re
import h5py
from datetime import datetime
import numpy as np
from tqdm import tqdm
from matplotlib import patches
from matplotlib import collections

logfiles = [f for f in os.listdir() if 'main_stage2.log' in f]

logs = []
for name in tqdm(logfiles):
    with open(name, 'r') as F:
        logs += F.readlines()


timepoints = {}
for line in tqdm(logs):
    ranknode, time, content = line.split(" | ")
    rank, node = ranknode.split(":")
    function, message = content.split(" >> ")
    time = datetime.fromisoformat( time.replace(',', '.') )
    try:
        timepoints[rank].append({'time':time.timestamp(), 'msg':message})
    except KeyError:
        timepoints[rank] = []
        timepoints[rank].append({'time':time.timestamp(), 'msg':message})

events = {}
for rank in tqdm(timepoints):
    iteration_start = []
    update_Fcell_start = []
    update_Fcell_end = []
    add_diffBragg_spots_start = []
    add_diffBragg_spots_end = []
    shots = {}
    MPI_barrier_starts = []
    MPI_barrier_end = []
    MPI_aggregation_starts = []
    MPI_aggregation_end = []
    for tp in timepoints[rank]:
        time = tp['time']
        msg = tp['msg']
        if 'BEGIN FUNC GRAD ; iteration' in msg:
            iteration_start.append(time)
        elif 'start update Fcell' in msg:
            update_Fcell_start.append(time)
        elif 'done update Fcell' in msg:
            update_Fcell_end.append(time)
        elif 'run diffBragg for shot' in msg:
            shot_id = int(msg.split('shot')[1])
            if shot_id in shots.keys():
                shots[shot_id]['start'].append(time)
            else:
                shots[shot_id] = {'start':[time], 'end':[]}
        elif 'finished diffBragg for shot' in msg:
            shot_id = int(msg.split('shot')[1])
            if shot_id in shots.keys():
                shots[shot_id]['end'].append(time)
            else:
                shots[shot_id] = {'start':[], 'end':[time]}
        elif 'Time rank worked on shots' in msg:
            shottime = float( msg.split('=')[1] )
            add_diffBragg_spots_start.append(time-shottime)
            add_diffBragg_spots_end.append(time)
            MPI_barrier_starts.append(time)
        elif 'MPI aggregation of func and grad' in msg:
            MPI_barrier_end.append(time)
            MPI_aggregation_starts.append(time)
        elif 'Time for MPIaggregation' in msg:
            MPI_aggregation_end.append(time)
    events[rank] = {'iteration':iteration_start,
                    'Fcell_start':update_Fcell_start,
                    'Fcell_end':update_Fcell_end,
                    'diffBragg_start':add_diffBragg_spots_start,
                    'diffBragg_end':add_diffBragg_spots_end,
                    'shots':shots,
                    'barrier_start':MPI_barrier_starts,
                    'barrier_end':MPI_barrier_end,
                    'aggregation_start':MPI_aggregation_starts,
                    'aggregation_end':MPI_aggregation_end}

with h5py.File('timestamps.h5', 'w') as F:
    for rank in tqdm(timepoints):
        timetable = events[rank]
        grp_rank = F.create_group(rank)
        grp_rank.create_dataset('iteration_start', data=np.array(timetable['iteration']))
        grp_rank.create_dataset('update_Fcell/start', data=np.array(timetable['Fcell_start']))
        grp_rank.create_dataset('update_Fcell/end', data=np.array(timetable['Fcell_end']))
        grp_rank.create_dataset('add_diffBragg_spots/start', data=np.array(timetable['diffBragg_start']))
        grp_rank.create_dataset('add_diffBragg_spots/end', data=np.array(timetable['diffBragg_end']))
        grp_rank.create_dataset('MPI_barrier/start', data=np.array(timetable['barrier_start']))
        grp_rank.create_dataset('MPI_barrier/end', data=np.array(timetable['barrier_end']))
        grp_rank.create_dataset('MPI_aggregation/start', data=np.array(timetable['aggregation_start']))
        grp_rank.create_dataset('MPI_aggregation/end', data=np.array(timetable['aggregation_end']))
        for shot_id in timetable['shots']:
            shot_start = timetable['shots'][shot_id]['start']
            shot_end = timetable['shots'][shot_id]['end']
            grp_rank.create_dataset('shots/%d/start'%shot_id, data=np.array(shot_start))
            grp_rank.create_dataset('shots/%d/end'%shot_id, data=np.array(shot_end))


sys.exit()

def get_color(r,g,b):
    return r/255.,g/255.,b/255.

def get_patches(tvals_start, tvals_stop, rank, tcolors):
    patch_list = []
    for i_t, (t1,t2) in enumerate(zip(tvals_start, tvals_stop)):
        xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
        patch = patches.Polygon(xy=xy, color=tcolors[i_t], closed=True) #, ec='k', lw=0.5)
        patch_list.append(patch)
    C = collections.PatchCollection(patch_list, match_original=True)
    return C

COLOR_green = get_color(105,183,100)
COLOR_yellow = get_color(255,221,113)
COLOR_red = get_color(242,108,100)
COLOR_blue = get_color(113,221,255)

ax = plt.gca()
num_ranks = len(timepoints)
max_t = max([events[r]['aggregation_end'][-1] for r in events])
for rank in tqdm(range(num_ranks)):
    timetable = events['RANK'+str(rank)]

    # plot Fcell update in yellow
    tstarts = timetable['Fcell_start']
    tstops = timetable['Fcell_end']
    tcolors = [COLOR_yellow]*len(timetable['Fcell_start'])

    # plot add_diffBragg_spots in red
    tstarts += timetable['diffBragg_start']
    tstops += timetable['diffBragg_end']
    tcolors += [COLOR_red]*len(timetable['diffBragg_start'])

    #plot MPI barriers in blue
    tstarts += timetable['barrier_start']
    tstops += timetable['barrier_end']
    tcolors += [COLOR_blue]*len(timetable['barrier_start'])

    #plot MPI barriers in green
    tstarts += timetable['aggregation_start']
    tstops += timetable['aggregation_end']
    tcolors += [COLOR_green]*len(timetable['aggregation_start'])

    coll = get_patches(tstarts, tstops, rank=rank, tcolors=tcolors)
    ax.add_collection(coll)
    ax.tick_params(labelsize=8, pad=1)
    ax.plot( timings['iteration'], [rank]*len(timings['iteration']), '*', mew=0, ms=4, color='k')

for r in range(num_ranks+1):
    ax.plot([0, max_t], [r-0.5, r-0.5], color='k', lw=0.33)

plt.ylim(-.5,num_ranks-0.5)
ax.set_xlabel("time after first iteration (seconds)", fontsize=12)
ax.set_ylabel("rank number", fontsize=12)
plt.xlim(0,max_t)
plt.gcf().set_size_inches((6,3))



patch_legend = [
    patches.Patch(color=COLOR_yellow, ec='k', lw=0.5,label="update_Fcell"),
    patches.Patch(color=COLOR_red, ec='k', lw=0.5, label='add_diffBragg_spots'),
    patches.Patch(color=COLOR_blue, ec='k', lw=0.5, label='MPI barrier'),
    patches.Patch(color=COLOR_green, ec='k', lw=0.5, label='MPI aggregation'),
    patches.Patch(color='w', ec='k', lw=0.5, label='untracked'),
    plt.plot([], [], ls="", marker='*', color='k',ms=5,mec=None,  label="Iteration begin")[0]
]

leg = ax.legend(handles=patch_legend, markerscale=1,
                 bbox_to_anchor=(.99,0.5),
                 prop={'size':7.5},
                 loc="center left")

fr = leg.get_frame()
fr.set_alpha(1)
fr.set_facecolor('w')
fr.set_edgecolor('w')

plt.subplots_adjust(right=0.8, bottom=0.15)

plt.show()
