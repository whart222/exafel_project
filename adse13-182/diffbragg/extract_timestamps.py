import os
import sys
import h5py
from datetime import datetime
import numpy as np
from tqdm import tqdm

logfiles = [f for f in os.listdir() if 'main_stage2.log' in f]

if len(logfiles)==0:
    print("Found no logfiles! Exiting.")
    sys.exit()

logs = []
for name in tqdm(logfiles):
    with open(name, 'r') as F:
        logs += F.readlines()


timepoints = {}
nodes = {}
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
        nodes[rank] = node


def get_template():
    return {'start':[], 'end':[]}

events = {}
for rank in tqdm(timepoints):
    startup = get_template()
    refinement = get_template()
    prep_dataframe = get_template()
    from_json_file = get_template()
    GatherFromExperiment = get_template()
    MPI_barrier_startup = get_template()
    gather_Hi_info = get_template()
    setup = get_template()
    iterations = get_template()
    update_Fcell = get_template()
    add_diffBragg_spots = get_template()
    shots = {}
    MPI_barrier_iterations = get_template()
    MPI_aggregation = get_template()

    for tp in timepoints[rank]:
        time = tp['time']
        msg = tp['msg']

        # A kingdom for a switch/case in python!
        if "EVENT: read input pickle" in msg:
            startup['start'].append(time)
        elif "_launch done run setup" in msg:
            startup['end'].append(time)

        elif "_launcher runno setup" in msg:
            refinement['start'].append(time)
        elif "_launcher done runno setup" in msg:
            refinement['end'].append(time)

        elif "EVENT: BEGIN prep dataframe" in msg:
            prep_dataframe['start'].append(time)
        elif "EVENT: DONE prep dataframe" in msg:
            prep_dataframe['end'].append(time)

        elif "EVENT: BEGIN loading experiment list" in msg:
            from_json_file['start'].append(time)
        elif "EVENT: DONE loading experiment list" in msg:
            from_json_file['end'].append(time)

        elif "EVENT: LOADING ROI DATA" in msg:
            GatherFromExperiment['start'].append(time)
        elif "EVENT: DONE LOADING ROI" in msg:
            GatherFromExperiment['end'].append(time)

        elif "DONE LOADING DATA; ENTER BARRIER" in msg:
            MPI_barrier_startup['start'].append(time)
        elif "DONE LOADING DATA; EXIT BARRIER" in msg:
            MPI_barrier_startup['end'].append(time)

        elif "EVENT: Gathering global HKL information" in msg:
            gather_Hi_info['start'].append(time)
        elif "EVENT: FINISHED gather global HKL information" in msg:
            gather_Hi_info['end'].append(time)

        elif "Setup begins!" in msg:
            setup['start'].append(time)
        elif "Setup ends!" in msg:
            setup['end'].append(time)

        elif "BEGIN FUNC GRAD ; iteration" in msg:
            iterations['start'].append(time)
        elif "DONE WITH FUNC GRAD" in msg:
            iterations['end'].append(time)

        elif "start update Fcell" in msg:
            update_Fcell['start'].append(time)
        elif "done update Fcell" in msg:
            update_Fcell['end'].append(time)

        elif "run diffBragg for shot" in msg:
            shot_id = int(msg.split('shot')[1])
            if shot_id not in shots.keys():
                shots[shot_id] = get_template()
            shots[shot_id]['start'].append(time)
        elif "finished diffBragg for shot" in msg:
            shot_id = int(msg.split('shot')[1])
            if shot_id not in shots.keys():
                shots[shot_id] = get_template()
            shots[shot_id]['end'].append(time)

        elif "Time rank worked on shots" in msg:
            shottime = float( msg.split('=')[1] )
            add_diffBragg_spots['start'].append(time-shottime)
            add_diffBragg_spots['end'].append(time)
            MPI_barrier_iterations['start'].append(time)
        elif 'MPI aggregation of func and grad' in msg:
            MPI_barrier_iterations['end'].append(time)
            MPI_aggregation['start'].append(time)
        elif "Time for MPIaggregation" in msg:
            MPI_aggregation['end'].append(time)

    events[rank] = {'startup': startup,
                    'refinement': refinement,
                    'prep_dataframe': prep_dataframe,
                    'from_json_file': from_json_file,
                    'GatherFromExperiment': GatherFromExperiment,
                    'MPI_barrier_startup': MPI_barrier_startup,
                    'gather_Hi_info': gather_Hi_info,
                    'setup': setup,
                    'iterations': iterations,
                    'update_Fcell': update_Fcell,
                    'add_diffBragg_spots': add_diffBragg_spots,
                    'shots': shots,
                    'MPI_barrier_iterations': MPI_barrier_iterations,
                    'MPI_aggregation': MPI_aggregation}

with h5py.File('timestamps.h5', 'w') as F:
    for rank in tqdm(timepoints):
        timetable = events[rank]
        grp_rank = F.create_group(rank)
        grp_rank.attrs['node'] = nodes[rank]

        for key in timetable:
            if key == 'shots':
                continue
            grp_rank.create_dataset(key+'/start', data=np.array(timetable[key]['start']))
            grp_rank.create_dataset(key+'/end', data=np.array(timetable[key]['end']))

        for shot_id in timetable['shots']:
            shot_start = timetable['shots'][shot_id]['start']
            shot_end = timetable['shots'][shot_id]['end']
            grp_rank.create_dataset('shots/%d/start'%shot_id, data=np.array(shot_start))
            grp_rank.create_dataset('shots/%d/end'%shot_id, data=np.array(shot_end))
