import os
#import re
import sys
import h5py
from datetime import datetime
import numpy as np
from tqdm import tqdm
from matplotlib import patches
from matplotlib import collections
from matplotlib import pyplot as plt

def get_color(r,g,b):
    return r/255.,g/255.,b/255.

#COLOR_green = get_color(105,183,100)
COLOR_azure = get_color(153, 201, 240)
COLOR_green = get_color(103, 191, 92)
COLOR_yellow = get_color(255,221,113)
COLOR_red = get_color(242,108,100)
COLOR_blue = get_color(23,190,207)
COLOR_orange = get_color(255,128,14)
COLOR_pink = get_color(247,182,210)
COLOR_gray = get_color(119,119,119)

def read_file(name):
    datapoints = {}
    with h5py.File(name, 'r') as F:
        for rank in tqdm(F['/']):
            times = {}
            if not 'RANK' in rank:
                continue
            for metric in F['/'+rank]:
                try:
                    times[metric+'_start'] = F['/'+rank+'/'+metric+'/start'][:]
                    times[metric+'_end'] = F['/'+rank+'/'+metric+'/end'][:]
                except:
                    pass
            times['shots'] = {}
            for shot in F['/'+rank+'/shots']:
                times['shots'][shot] = {}
                times['shots'][shot]['start'] = F['/'+rank+'/shots/'+shot+'/start'][:]
                times['shots'][shot]['end'] = F['/'+rank+'/shots/'+shot+'/end'][:]
            datapoints[rank] = times
    return datapoints

def adjust_starttime(datapoints):
    starttime = min( [datapoints[rank]['refinement_start'][0] for rank in datapoints] )
    for rank in datapoints:
        for metric in datapoints[rank]:
            if type(datapoints[rank][metric])==np.ndarray:
                datapoints[rank][metric] -= starttime
        for shot in datapoints[rank]['shots']:
            datapoints[rank]['shots'][shot]['start'] -= starttime
            datapoints[rank]['shots'][shot]['end'] -= starttime
    return starttime

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

def print_iteration_times(datapoints):
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

def add_total_time(timetable, metric):
    return sum(timetable[metric+'end']-timetable[metric+'start'])

def print_times(datapoints):
    num_ranks = len(datapoints)
    #t_zero = min([datapoints[r]['iteration'][0] for r in datapoints])
    #t_end = max([datapoints[r]['aggregation_end'][-1] for r in datapoints])
    #total_time = (t_end - t_zero) * num_ranks

    t_Fcell = []
    t_diffbragg = []
    t_barrier = []
    t_aggregation = []
    t_iterations = []

    for rank in tqdm(range(num_ranks)):
        timetable = datapoints['RANK'+str(rank)]

        t_Fcell.append( add_total_time(timetable, 'update_Fcell_') )
        t_diffbragg.append(0)
        for shot in timetable['shots']:
            t_diffbragg[-1] += add_total_time(timetable['shots'][shot], '')
        t_barrier.append( add_total_time(timetable, 'MPI_barrier_iterations_') )
        t_aggregation.append( add_total_time(timetable, 'MPI_aggregation_') )
        t_iterations.append( add_total_time(timetable, 'iterations_'))


    print("Time for refinement:", np.mean(t_iterations), 'sec')
    print("Time to update_Fcell:", np.mean(t_Fcell), 'sec')
    print("Time to add_diffBragg:", np.mean(t_diffbragg), 'sec')
    print("Time for MPI barrier:", np.mean(t_barrier), 'sec')
    print("Time for MPI aggregate:", np.mean(t_aggregation), 'sec')


class time_patches(object):

    def __init__(self):
        self.reset()

    def reset(self):
        self.tstarts = []
        self.tstops = []
        self.tcolors = []

    def add_patches(self, starts, stops, color):
        assert len(starts)==len(stops), 'ERROR! Number of starts=%d different from stops=%d'%(len(starts), len(stops))
        self.tstarts += list(starts)
        self.tstops += list(stops)
        self.tcolors += [color]*len(starts)

    def get_patches(self, rank):
        patch_list = []
        for i_t, (t1,t2) in enumerate(zip(self.tstarts, self.tstops)):
            xy = np.array([(t1,rank-0.5), (t1, rank+0.5), (t2, rank+0.5), (t2, rank-0.5)])
            #breakpoint()
            patch = patches.Polygon(xy=xy, color=self.tcolors[i_t], closed=True) #, ec='k', lw=0.5)
            patch_list.append(patch)
        C = collections.PatchCollection(patch_list, match_original=True)
        return C


def plot_startup(datapoints, limit_ranks=None):
    ax = plt.gca()
    if limit_ranks is None:
        num_ranks = len(datapoints)
    else:
        num_ranks = limit_ranks
    min_t = min([datapoints[r]['startup_start'][0] for r in datapoints])
    timepatch = time_patches()
    for rank in tqdm(range(num_ranks)):
        timetable = datapoints['RANK'+str(rank)]
        timepatch.reset()

        timepatch.add_patches(timetable['prep_dataframe_start'],
                          timetable['prep_dataframe_end'],
                          COLOR_azure)

        timepatch.add_patches(timetable['from_json_file_start'],
                          timetable['from_json_file_end'],
                          COLOR_gray)

        timepatch.add_patches(timetable['GatherFromExperiment_start'],
                          timetable['GatherFromExperiment_end'],
                          COLOR_green)

        timepatch.add_patches(timetable['MPI_barrier_startup_start'],
                          timetable['MPI_barrier_startup_end'],
                          COLOR_blue)

        timepatch.add_patches(timetable['gather_Hi_info_start'],
                          timetable['gather_Hi_info_end'],
                          COLOR_orange)

        timepatch.add_patches(timetable['setup_start'],
                          timetable['setup_end'],
                          COLOR_pink)

        coll = timepatch.get_patches(rank)
        ax.add_collection(coll)
    ax.tick_params(labelsize=8, pad=1)
    ax.plot( [min_t, min_t], [-.5,num_ranks-0.5], '--', lw=2, color='k')

    if limit_ranks is not None:
        for r in range(num_ranks+1):
            ax.plot([min_t, 0], [r-0.5, r-0.5], color='k', lw=0.33)

    plt.ylim(-.5,num_ranks-0.5)
    ax.set_xlabel("time before first iteration (seconds)", fontsize=12)
    ax.set_ylabel("rank number", fontsize=12)
    plt.xlim(1.1*min_t,0)
    plt.gcf().set_size_inches((12,6))

    patch_legend = [
        patches.Patch(color=COLOR_azure, ec='k', lw=0.5,label="prep_dataframe"),
        patches.Patch(color=COLOR_gray, ec='k', lw=0.5,label="from_json_file"),
        patches.Patch(color=COLOR_green, ec='k', lw=0.5, label='GatherFromExperiment'),
        patches.Patch(color=COLOR_blue, ec='k', lw=0.5, label='MPI barrier'),
        patches.Patch(color=COLOR_orange, ec='k', lw=0.5, label='gather_Hi_info'),
        patches.Patch(color=COLOR_pink, ec='k', lw=0.5, label='setup'),
        patches.Patch(color='w', ec='k', lw=0.5, label='untracked'),
        plt.plot([], [], ls="--", lw=2, color='k', label="program start")[0]
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


def plot_overview(datapoints, limit_ranks=None, limit_time=None, combine_shots=True):
    ax = plt.gca()
    if limit_ranks is None:
        num_ranks = len(datapoints)
    else:
        num_ranks = limit_ranks
    max_t = max([datapoints[r]['iterations_end'][-1] for r in datapoints])
    timepatch = time_patches()
    for rank in tqdm(range(num_ranks)):
        timetable = datapoints['RANK'+str(rank)]
        timepatch.reset()

        timepatch.add_patches(timetable['update_Fcell_start'],
                              timetable['update_Fcell_end'],
                              COLOR_yellow)

        if combine_shots:
            timepatch.add_patches(timetable['add_diffBragg_spots_start'],
                                  timetable['add_diffBragg_spots_end'],
                                  COLOR_red)
        else:
            for shot in timetable['shots']:
                timepatch.add_patches(timetable['shots'][shot]['start'],
                                      timetable['shots'][shot]['end'],
                                      COLOR_red)

        timepatch.add_patches(timetable['MPI_barrier_iterations_start'],
                              timetable['MPI_barrier_iterations_end'],
                              COLOR_blue)

        timepatch.add_patches(timetable['MPI_aggregation_start'],
                              timetable['MPI_aggregation_end'],
                              COLOR_green)

        coll = timepatch.get_patches(rank)
        ax.add_collection(coll)

        if limit_time is not None:
            ax.plot( list(timetable['iterations_start']), [rank]*len(timetable['iterations_start']), 'o', mew=0, ms=4, color='k')

    if limit_ranks is not None:
        for r in range(num_ranks+1):
            ax.plot([0, max_t], [r-0.5, r-0.5], color='k', lw=0.33)

    plt.ylim(-.5,num_ranks-0.5)
    ax.set_xlabel("time after first iteration (seconds)", fontsize=12)
    ax.set_ylabel("rank number", fontsize=12)
    if limit_time is None:
        plt.xlim(0,max_t)
    else:
        plt.xlim(0,limit_time)
    ax.tick_params(labelsize=8, pad=1)
    plt.gcf().set_size_inches((12,6))


    patch_legend = [
        patches.Patch(color=COLOR_yellow, ec='k', lw=0.5,label="update_Fcell"),
        patches.Patch(color=COLOR_red, ec='k', lw=0.5, label='add_diffBragg_spots'),
        patches.Patch(color=COLOR_blue, ec='k', lw=0.5, label='MPI barrier'),
        patches.Patch(color=COLOR_green, ec='k', lw=0.5, label='MPI aggregation'),
        patches.Patch(color='w', ec='k', lw=0.5, label='untracked'),
        plt.plot([], [], ls="", marker='o', color='k',ms=5,mec=None,  label="Iteration begin")[0]
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
    if len(sys.argv)>1:
        datapoints = read_file(sys.argv[1])
        adjust_starttime(datapoints)
        print_times(datapoints)
    #plot_overview(datapoints)
