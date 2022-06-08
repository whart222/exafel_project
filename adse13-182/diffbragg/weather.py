import os
#import re
import h5py
from datetime import datetime
import numpy as np
from tqdm import tqdm
from matplotlib import patches
from matplotlib import collections
from matplotlib import pyplot as plt

def get_color(r,g,b):
    return r/255.,g/255.,b/255.

COLOR_green = get_color(105,183,100)
COLOR_yellow = get_color(255,221,113)
COLOR_red = get_color(242,108,100)
COLOR_blue = get_color(113,221,255)

def get_patches(tvals_start, tvals_stop, rank, tcolors):
    patch_list = []
    for i_t, (t1,t2) in enumerate(zip(tvals_start, tvals_stop)):
        xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
        patch = patches.Polygon(xy=xy, color=tcolors[i_t], closed=True) #, ec='k', lw=0.5)
        patch_list.append(patch)
    C = collections.PatchCollection(patch_list, match_original=True)
    return C

def read_file(name):
    datapoints = {}
    with h5py.File(name, 'r') as F:
        for rank in tqdm(F['/']):
            times = {}
            if not 'RANK' in rank:
                continue
            times['iteration'] = F['/'+rank+'/iteration_start'][:]
            times['Fcell_start'] = F['/'+rank+'/update_Fcell/start'][:]
            times['Fcell_end'] = F['/'+rank+'/update_Fcell/end'][:]
            times['diffBragg_start'] = F['/'+rank+'/add_diffBragg_spots/start'][:]
            times['diffBragg_end'] = F['/'+rank+'/add_diffBragg_spots/end'][:]
            times['barrier_start'] = F['/'+rank+'/MPI_barrier/start'][:]
            times['barrier_end'] = F['/'+rank+'/MPI_barrier/end'][:]
            times['aggregation_start'] = F['/'+rank+'/MPI_aggregation/start'][:]
            times['aggregation_end'] = F['/'+rank+'/MPI_aggregation/end'][:]
            datapoints[rank] = times
    return datapoints


class online_stats(object):
    def __init__(self, unit=""):
        self.count = 0
        self.mean = 0
        self.M2 = 0
        self.unit = ""

    def add(self, number):
        # Use temporaries to avoid propagating errors
        new_count = self.count + 1
        delta = number - self.mean
        new_mean = self.mean + delta / new_count
        new_M2 = M2 + delta*(number - new_mean)

        self.count = count
        self.mean = new_mean
        self.M2 = new_M2

    def add_multi(self, count, mean, var):
        new_count = self.count + count
        delta = mean - self.mean
        new_mean = (mean*count + self.mean*self.count) / new_count
        new_M2 = count*var + self.M2 + delta**2 * self.count*count / new_count

        self.count = new_count
        self.mean = new_mean
        self.M2 = new_M2

    def add_ndarray(self, data):
        self.add_multi(data.size, data.mean(), data.var())

    def get_mean(self):
        return self.mean

    def get_sum(self):
        return self.mean * self.count

    def get_var(self):
        return self.M2 / self.count

    def get_std(self):
        return self.get_var()**0.5

    def __str__(self):
        return f"{self.mean}{self.unit} +/- {self.get_std()}{self.unit}"

def print_times(datapoints):
    num_ranks = len(datapoints)
    t_zero = min([datapoints[r]['iteration'][0] for r in datapoints])
    t_end = max([datapoints[r]['aggregation_end'][-1] for r in datapoints])
    total_time = (t_end - t_zero) * num_ranks

    t_Fcell = online_stats("sec")
    t_diffbragg = online_stats("sec")
    t_barrier = online_stats("sec")
    t_aggregation = online_stats("sec")

    for rank in tqdm(range(num_ranks)):
        timetable = datapoints['RANK'+str(rank)]

        t_Fcell.add_ndarray(timetable['Fcell_end'] - timetable['Fcell_start'])
        t_diffbragg.add_ndarray(timetable['diffBragg_end'] - timetable['diffBragg_start'])
        t_barrier.add_ndarray(timetable['barrier_end'] - timetable['barrier_start'])
        t_aggregation.add_ndarray(timetable['aggregation_end'] - timetable['aggregation_start'])

    print("Time to update_Fcell:", t_Fcell)
    print("Time to add_diffBragg:", t_diffbragg)
    print("Time for MPI barrier:", t_barrier)
    print("Time for MPI aggregate:", t_aggregation)


def plot_overview(datapoints):
    ax = plt.gca()
    num_ranks = len(datapoints)
    t_zero = datapoints['RANK0']['iteration'][0]
    max_t = max([datapoints[r]['aggregation_end'][-1] for r in datapoints]) - t_zero
    for rank in tqdm(range(num_ranks)):
        timetable = datapoints['RANK'+str(rank)]

        # set T0 to first iteration
        for listname in timetable:
            timetable[listname] -= t_zero

        # plot Fcell update in yellow
        tstarts = list(timetable['Fcell_start'])
        tstops = list(timetable['Fcell_end'])
        tcolors = [COLOR_yellow]*len(timetable['Fcell_start'])

        # plot add_diffBragg_spots in red
        tstarts += list(timetable['diffBragg_start'])
        tstops += list(timetable['diffBragg_end'])
        tcolors += [COLOR_red]*len(timetable['diffBragg_start'])

        #plot MPI barriers in blue
        tstarts += list(timetable['barrier_start'])
        tstops += list(timetable['barrier_end'])
        tcolors += [COLOR_blue]*len(timetable['barrier_start'])

        #plot MPI barriers in green
        tstarts += list(timetable['aggregation_start'])
        tstops += list(timetable['aggregation_end'])
        tcolors += [COLOR_green]*len(timetable['aggregation_start'])

        coll = get_patches(tstarts, tstops, rank=rank, tcolors=tcolors)
        ax.add_collection(coll)
        ax.tick_params(labelsize=8, pad=1)
        #ax.plot( list(timetable['iteration']), [rank]*len(timetable['iteration']), '*', mew=0, ms=4, color='k')

    #for r in range(num_ranks+1):
    #    ax.plot([0, max_t], [r-0.5, r-0.5], color='k', lw=0.33)

    plt.ylim(-.5,num_ranks-0.5)
    ax.set_xlabel("time after first iteration (seconds)", fontsize=12)
    ax.set_ylabel("rank number", fontsize=12)
    plt.xlim(0,max_t)
    plt.gcf().set_size_inches((12,6))



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

    #plt.subplots_adjust(right=0.8, bottom=0.15)

    plt.show()




if __name__=='__main__':
    datapoints = read_file('timestamps.h5')
    plot_overview(datapoints)
