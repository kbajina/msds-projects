"""
Northwestern University
MSDS 460 - Decision Analytics

Home Building Project Simulation
"""

from sympy.utilities.iterables import multiset_permutations
import numpy as np
import pandas as pd
pd.set_option('display.expand_frame_repr', False)
from collections import defaultdict

# =============================================================================
# Identify Valid Project Traces (Sequence of Activities)
# =============================================================================

# identify all project activities
activities = ['dig','bricks','tiles', 'lay', 'walls', 
              'wiring', 'roofing', 'painting']

# list of lists to track all possible project traces
all_traces = []

# lists of lists to track all valid traces (given problem constraints)
valid_traces = []

# iterate through all permutations of activities
# from these permutations we will drop all traces (scenarios) that do not
# meet our sequence constraints (as defined in the problem).
for trace in multiset_permutations(activities):
    all_traces.append(trace)
    
    # identify order of activities in the trace (using index value)
    dig_index = trace.index('dig') 
    bricks_index = trace.index('bricks') 
    tiles_index = trace.index('tiles') 
    lay_index = trace.index('lay')
    walls_index = trace.index('walls') 
    wiring_index = trace.index('wiring')
    roofing_index = trace.index('roofing')
    painting_index = trace.index('painting')
    
    # only capture traces that meet sequence constraints
    # project trace (sequence) can begin with Dig, Brick or Tiles
    # Dig must occur BEFORE Lay
    # Bricks must occur BEFORE walls
    # Walls must occur BEFORE wiring
    # Tiles must occur BEFORE roofing
    # Walls must occur BEFORE roofing
    # Wiring must occur BEFORE painting
    # Roofing must occur BEFORE painting
    if (dig_index == 0 or bricks_index == 0 or tiles_index == 0)\
        and (dig_index < lay_index)\
        and (lay_index < walls_index)\
        and (bricks_index < walls_index)\
        and (walls_index < wiring_index)\
        and (tiles_index < roofing_index)\
        and (walls_index < roofing_index)\
        and (wiring_index < painting_index)\
        and (roofing_index < painting_index):
            valid_traces.append(trace)

all_traces

# =============================================================================
# Define Service Time Distributions
# =============================================================================

# list modal times (in order of "activities" list)
times = [4, # dig
         7, # bricks
         12, # tiles
         2, # lay
         10, # walls
         3, # wiring
         5, # roofing
         4] # painting

# dictionary of all model times associated with an activity
activity_times = dict(zip(activities, times))

# dictionary of all minimum times for each activity
# min time = 1/2 * modal times
activities_min = defaultdict(int)
for act,time in activity_times.items():
    activities_min[act] = time * 1/2
activities_min

# dictionary of all maximum times for each activity
# max time = 2 * modal times
activities_max = defaultdict(int)
for act,time in activity_times.items():
    activities_max[act] = time * 2
activities_max


# use random triangular distribution for event/activity times
def random_triangle_sample(min:float, mode:float, max:float):
    """

    Parameters
    ----------
    min : float
        Minimum value in the distribution for sampling. The default is float.
    mode : float
        Mode (most frequently occuring) value in the distribution for sampling.
        The default is float.
    max : float
        Maximum value in the distribution for sampling. The default is float.

    Returns
    -------
    random_sample: float.

    """
    random_sample = np.random.triangular(min, mode, max)
    
    return random_sample

# NOTE: time for an event/activity should be reduced in proportion to the 
# number of workers assigned to it.
def total_activity_time(t:float, w:int):
    """
    
    Parameters
    ----------
    t : float
        Total time expected for an activity, assuming ONE worker.
    w : int
        Total workers available for the particular activity.

    Returns
    -------
    total_time: float.

    """
    # assumption is that all activity times are based on TWO workers
    # divide expected time for an activity by the number of available workers
    # if less than 2 workers are used, time will increase (proportionally)
    # if more than 2 workers are used, time will decrease (proportionally)
    total_time = t * (2/ w)
    
    return total_time


# =============================================================================
# Worker Resource Availability
# =============================================================================

# determine worker requirements by activity
worker_req = [1, # dig
              1, # bricks
              1, # tiles
              1, # lay
              1, # walls
              2, # wiring
              2, # roofing
              1] # painting

# create dictionary of minimum workers needed for each activity
activity_min_workers = dict(zip(activities, worker_req))

# =============================================================================
# Monte Carlo Simulation of Home Building Project
# Output: Project Case-by-Case Event Log
# =============================================================================

    # Dig must occur BEFORE Lay
    # Bricks must occur BEFORE walls
    # Walls must occur BEFORE wiring
    # Tiles must occur BEFORE roofing
    # Walls must occur BEFORE roofing
    # Wiring must occur BEFORE painting
    # Roofing must occur BEFORE painting

# identify which activities can be done concurrently
concurrent_act = {'dig': ['bricks', 'tiles'],
                  'bricks': ['dig', 'tiles'],
                  'tiles' : ['dig', 'bricks'],
                  'lay': ['bricks', 'tiles'],
                  'walls':['tiles', 'roofing'],
                  'wiring':['tiles', 'roofing'],
                  'roofing':['wiring']
    }

# identify total iterations for simulation
iterations = 10_000

# create blank dataframe to append minimum trace times from each iteration
min_tracker = pd.DataFrame(columns=["sens_w_count", "iteration",
                                    "caseid", "time"])


# min trace dict for each iteration
min_trace_ids = defaultdict(list)

# Sensitivity Analysis: determine potential range of worker count
min_workers = 2
max_workers = 10
workers_ranges = range(min_workers,max_workers+1)

# Full set of "iterations" will be run for every range of workers.
# Each iteration will evaluate the potential completion time for all 33 valid
# project completion traces.
# This analysis will be used to determine how the optimal trace might change 
# given more workers available to complete project activities.

for workers in workers_ranges:
    print('\n')
    print(f'Total Workers: {workers}')
    print("--------------------------")
    
    for iteration in range(iterations):
        # track output from each trace (case id)
        trace_dict = defaultdict(list)
            
        # iterate through all valid traces (sequences) of activities
        for trace in valid_traces:
            
            # track workers used/ available during each trace
            worker_count = workers
            
            # track timestamp for beginning and end of each activity
            time_stamp = 0
        
            # Append starting point for each trace
            trace_dict['caseid'].append(valid_traces.index(trace))
            trace_dict['time'].append(0)
            trace_dict['activity'].append('start')
            trace_dict['act_type'].append('start')
            trace_dict['workers_used'].append(None)
            trace_dict['time_worked'].append(None)
            trace_dict['workers_avail'].append(worker_count)
        
            
            # create variable to monitor if concurrent activities were performed
            # in such an event, details of an activity were defined in the previous
            # iteration of looping through an activity in the trace.
            # this variable will indicate if we must skip the ith activity
            skip_act = False
            skip_count = 0
        
            
            # iterate through every index in trace
            # current index is ith activity in trace list
            # next index is the following activity in trace list
            for i in range(len(trace)):
                # check if skip has already taken place (i.e. skip_count==2)
                if skip_act and skip_count==1:
                    skip_count += 1
                    continue # skip to next iteration of the loop
                    
                # reset flag for skipping an activity in the loop
                skip_act = False
                skip_count = 0
                    
                if i < len(trace)-1:
                    current_act = trace[i]
                    next_act = trace[i+1]
                else:
                    current_act = trace[i]
                    next_act = None
                
                # =====================================================================
                # Beginning of Current Activity
                # =====================================================================
                
                # track start time of activity
                c_act_start = time_stamp
                
                # track number of workers available before activities started
                # will not be a modified counter like "worker_count"
                start_avail_workers = worker_count
                
                # Append current activity beginning stats to trace dictionary
                trace_dict['caseid'].append(valid_traces.index(trace))
                trace_dict['time'].append(c_act_start)
                trace_dict['activity'].append('beg_' + current_act)
                trace_dict['act_type'].append('primary')
                trace_dict['workers_used'].append(None)
                trace_dict['time_worked'].append(None)
                trace_dict['workers_avail'].append(start_avail_workers)
        
                # =====================================================================
                # Assess Requirements for Current Activity
                # =====================================================================
                
                # Determine projected activity time 
                # activity time distribution
                c_min = activities_min[current_act]
                c_mode = activity_times[current_act]
                c_max = activities_max[current_act]
                
                # sample from triangle distribution
                c_proj_time = random_triangle_sample(c_min, c_mode, c_max)
                
                # determine how many workers are needed
                c_req_workers = activity_min_workers[current_act]
        
                # start with required workers; updated further down if more available
                c_workers_used = c_req_workers
                
                # determine final time based on workers used
                c_total_time = total_activity_time(c_proj_time, c_req_workers)
                
                # track end time of activity
                c_end_time = time_stamp + c_total_time
                
                # reduce the available worker count
                worker_count += -(c_req_workers)
                
                # =====================================================================
                # Concurrent Work
                # =====================================================================
                # determine if the next activity in the trace can be done concurrently
                # if it cannot, can we use remaining workers to speed up the current
                # activity completion time?
                
                # only continue if we are NOT on the last activity in the trace, AND
                # if next activity is a valid concurrent activity
                if i < len(trace)-1 and (next_act in concurrent_act[current_act]):
                    # determine how many workers are needed for NEXT activity
                    n_req_workers = activity_min_workers[next_act]
                    
                    # determine if required number of workers are available
                    if n_req_workers <= worker_count:
                        
                        # =============================================================
                        # Beginning of Next (Concurrent) Activity
                        # =============================================================
                        
                        # Append next activity beginning stats to trace dictionary
                        trace_dict['caseid'].append(valid_traces.index(trace))
                        trace_dict['time'].append(c_act_start) # starts with current act.
                        trace_dict['activity'].append('beg_' + next_act)
                        trace_dict['act_type'].append('concurrent')
                        trace_dict['workers_used'].append(None)
                        trace_dict['time_worked'].append(None)
                        trace_dict['workers_avail'].append(start_avail_workers)
        
                    
                        # =============================================================
                        # Ending of Next (Concurrent) Activity
                        # =============================================================
                        
                        # Determine projected time for NEXT activity 
                        # activity time distribution
                        n_min = activities_min[next_act]
                        n_mode = activity_times[next_act]
                        n_max = activities_max[next_act]
                    
                        # sample from triangle distribution
                        n_proj_time = random_triangle_sample(n_min, n_mode, n_max)
                        
                        # if projected time is longer for next activity than current
                        # activity, put remaining resources towards next activity
                        if n_proj_time > c_proj_time:
                            n_req_workers = worker_count
                        
                        # determine final time based on workers used
                        n_total_time = total_activity_time(n_proj_time, n_req_workers)
                        
                        # track end time of next (concurrent) activity
                        n_end_time = c_act_start + n_total_time
                        
                        # reduce the available worker count
                        worker_count += -(n_req_workers)
                        
                        # remaining workers
                        
                        # flag that the next activity should be skipped in the loop
                        # since the work will be done concurrently with current activity
                        skip_act = True
                        skip_count += 1
                        
                    # else:
                    #     # don't skip next activity in the loop
                    #     skip_act = False
                        
                # =====================================================================
                # Updated Completion Info for Current Activity
                # =====================================================================
                # utilize remaining workers (whether concurrent activity was
                # completed or not) to improve speed of completion for current
                # activity in trace
                if worker_count > 0:
                    # update current activity time based on required workers
                    # plus number of remaining workers
                    c_total_time = total_activity_time(c_proj_time, 
                                                       c_req_workers + worker_count)
                    
                    # updated end time of activity
                    c_end_time = c_act_start + c_total_time
                    
                    # track total workers used for current activity
                    c_workers_used = c_req_workers + worker_count
                    
                    # reduce the available worker count to zero now that all are used
                    worker_count = 0
                    
                    
                # =====================================================================
                # Record Current and Next (Concurrent) Activity Completion Info
                # =====================================================================
                
                # if concurrent work was completed, determine which should be 
                # recorded first in the trace dictionary
                if skip_act:
                    # record the shorter of the two activities first
                    if c_end_time > n_end_time:
                        
                        # Append next (concurrent) activity ending stats to dict
                        trace_dict['caseid'].append(valid_traces.index(trace))
                        trace_dict['time'].append(n_end_time)
                        trace_dict['activity'].append('end_' + next_act)
                        trace_dict['act_type'].append('concurrent')
                        trace_dict['workers_used'].append(n_req_workers)
                        trace_dict['time_worked'].append(n_total_time)
                        trace_dict['workers_avail'].append(start_avail_workers
                                                           - n_req_workers)
        
                        
                        # Append current activity ending stats to dict
                        trace_dict['caseid'].append(valid_traces.index(trace))
                        trace_dict['time'].append(c_end_time)
                        trace_dict['activity'].append('end_' + current_act)
                        trace_dict['act_type'].append('primary')
                        trace_dict['workers_used'].append(c_workers_used)
                        trace_dict['time_worked'].append(c_total_time)
                        trace_dict['workers_avail'].append(worker_count)
        
                        
                        # reset worker count before next loop (activity)
                        worker_count = (worker_count
                                        + n_req_workers
                                        + c_workers_used)
                        
                        # update final timestamp
                        time_stamp = c_end_time
                        
                    # swith order of info recorded if c_end_time < n_end_time
                    else:
                        # Append current activity ending stats to dict
                        trace_dict['caseid'].append(valid_traces.index(trace))
                        trace_dict['time'].append(c_end_time)
                        trace_dict['activity'].append('end_' + current_act)
                        trace_dict['act_type'].append('primary')
                        trace_dict['workers_used'].append(c_workers_used)
                        trace_dict['time_worked'].append(c_total_time)
                        trace_dict['workers_avail'].append(start_avail_workers
                                                           - c_workers_used)
        
                        
                        # Append next (concurrent) activity ending stats to dict
                        trace_dict['caseid'].append(valid_traces.index(trace))
                        trace_dict['time'].append(n_end_time)
                        trace_dict['activity'].append('end_' + next_act)
                        trace_dict['act_type'].append('concurrent')
                        trace_dict['workers_used'].append(n_req_workers)
                        trace_dict['time_worked'].append(n_total_time)
                        trace_dict['workers_avail'].append(worker_count)
                        
                        # reset worker count before next loop (activity)
                        worker_count = (worker_count
                                        + n_req_workers
                                        + c_workers_used)
                        
                        # update final timestamp
                        time_stamp = n_end_time
    
                else:
                    # In the event of no concurrent work being completed, only
                    # append current activity ending stats to trace dictionary
                    trace_dict['caseid'].append(valid_traces.index(trace))
                    trace_dict['time'].append(c_end_time)
                    trace_dict['activity'].append('end_' + current_act)
                    trace_dict['act_type'].append('primary')
                    trace_dict['workers_used'].append(c_workers_used)
                    trace_dict['time_worked'].append(c_total_time)
                    trace_dict['workers_avail'].append(start_avail_workers
                                                       - c_workers_used)
                    
                    # reset worker count before next loop (activity)
                    worker_count = start_avail_workers
                
                    # update final timestamp
                    time_stamp = c_end_time
        
        # dataframe summarizing outcome of all traces
        traces_df = pd.DataFrame(trace_dict)
        
        # identify final times for each case trace
        case_times = traces_df.groupby(["caseid"])["time"].max()
        
        case_times = pd.DataFrame(case_times).reset_index(level=0)
        case_times['iteration'] = pd.Series([iteration +1]*33)
        case_times['sens_w_count'] = pd.Series([workers]*33)
        case_times = case_times[['sens_w_count', 'iteration',
                                 'caseid', 'time']]
        
        # append iteration results to main tracker
        min_tracker= pd.concat([min_tracker, case_times], axis=0)
        
        # determine trace ID with min total time
        min_trace = case_times.loc[case_times['time']==case_times['time'].min(),
                                   'caseid'].values[0]
        
        # append shortest path to tracker
        min_trace_ids['sens_w_count'].append(workers)
        min_trace_ids['iteration'].append(iteration +1)
        min_trace_ids['min_trace'].append(min_trace)
        
        if iteration % 10 == 0:
            print('', 
                  end=f'\rIteration {iteration+10}/{iterations} completed.')
        
        

# =============================================================================
# Determine the most frequently occuring traces, that have the lowest total
# activity time across all simulation iterations, for each level of 
# workers available.
# =============================================================================

# create pandas dataframe for min traces and their frequencies
min_trace_ids_df = pd.DataFrame(min_trace_ids)

# store summary of sensitivity analysis for most frequently occuring min trace
sens_min_tracker_df = pd.DataFrame(columns=["sens_w_count", "min_trace", "iteration"])

for workers in workers_ranges:
    # determine valid trace with the maximal occurences of having min time
    worker_min_df = min_trace_ids_df.loc[min_trace_ids_df['sens_w_count']==workers,]
    
    optimal_trace = worker_min_df.groupby(["sens_w_count","min_trace"])['iteration'].count().reset_index()
    
    optimal_trace = optimal_trace \
                        .loc[optimal_trace['iteration']==optimal_trace['iteration'].max(),]  \
                        .reset_index(drop=True)
                 
    # append iteration results to sensitivity tracker
    sens_min_tracker_df = pd.concat([sens_min_tracker_df, 
                                        optimal_trace], axis=0)

sens_min_tracker_df.rename(columns={"sens_w_count": "workers", 
                                    "min_trace": "trace_id",
                                    "iteration": "iter_count"},
                           inplace=True)



# =============================================================================
# Determine the trace with minimum mean time across all simulation iterations,
# for each level of workers available.
# =============================================================================

# store summary of sensitivity analysis for trace with smallest mean total time
sens_mean_tracker_df = pd.DataFrame(columns=["sens_w_count", "caseid", "time"])


for workers in workers_ranges:
    # determine mean times across all simulations for each caseid (trace)
    worker_mean_df = min_tracker.loc[min_tracker['sens_w_count']==workers,]
    min_mean_traces = worker_mean_df.groupby(['sens_w_count','caseid'])['time'].mean().reset_index()
    min_mean_traces = min_mean_traces \
                            .loc[min_mean_traces['time']==min_mean_traces['time'].min(),]  \
                            .reset_index(drop=True)

    # append iteration results to sensitivity tracker
    sens_mean_tracker_df = pd.concat([sens_mean_tracker_df, 
                                        min_mean_traces], axis=0)

sens_mean_tracker_df.rename(columns={"sens_w_count": "workers", 
                                    "caseid": "trace_id",
                                    "time": "mean_time"},
                           inplace=True)

# =============================================================================
# View outcomes to determine trace (sequence of activities, including 
# concurrent work, that will most likely yeild the best outcome for the 
# home building project.
# =============================================================================

sens_min_tracker_df

sens_mean_tracker_df
