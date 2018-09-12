#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd
from IPython import embed as shell

from sim_tools import DDM_traces_get, one_accumulater_traces_apply_bounds
from sim_tools import summary_plot, conditional_response_plot

sns.set(style='ticks', font='Arial', font_scale=1, rc={
    'axes.linewidth': 0.25, 
    'axes.labelsize': 7, 
    'axes.titlesize': 7, 
    'xtick.labelsize': 6, 
    'ytick.labelsize': 6, 
    'legend.fontsize': 6, 
    'xtick.major.width': 0.25, 
    'ytick.major.width': 0.25,
    'text.color': 'Black',
    'axes.labelcolor':'Black',
    'xtick.color':'Black',
    'ytick.color':'Black',} )
sns.plotting_context()

# data_folder = '/home/degee/research/model_simulations/ddm_data/'
# fig_folder = '/home/degee/research/model_simulations/ddm_figs/'
data_folder = '/Users/janwillem/Desktop/simulations/ddm_data/'
fig_folder = '/Users/janwillem/Desktop/simulations/ddm_figs/'

t = 5
dt = 0.01
timesteps = int(t/dt)
print(timesteps)
    
def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    for stim in [1,0]:
        x1 = DDM_traces_get( v=params['v'],
                                z=params['z'],
                                dc=params['dc'],
                                dc_slope=params['dc_slope'],
                                stim=stim,
                                nr_trials=params['nr_trials'],
                                dt=dt,)

        rt_dum, response_dum = one_accumulater_traces_apply_bounds(x1, 
                                a=params['a'],
                                b0_collapse=params['b0_collapse'],
                                b1_collapse=params['b1_collapse'],
                                dt=dt,)
        rt.append(rt_dum)
        response.append(response_dum)
        stimulus.append(np.ones(params['nr_trials']) * stim)

    df = pd.DataFrame()
    df.loc[:,'rt'] = (np.concatenate(rt) * dt) + ndt
    df.loc[:,'response'] = np.concatenate(response)
    df.loc[:,'stimulus'] = np.concatenate(stimulus)
    df.loc[:,'correct'] = np.array(np.concatenate(stimulus) == np.concatenate(response), dtype=int)
    df.loc[:,'subj_idx'] = params['subj_idx']
    df.to_csv(os.path.join(data_folder, 'df_{}.csv'.format(params['subj_idx'])))

simulate = False
nr_trials = int(1e5) #100K

v = 0.1
a = 0.1
ndt = 0.1
dc = 0
dc_slope = 0
ou = 0

sArray = [
    
    # 0 DDM neutral
    {'subj_idx':0, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},

    # 1 DDM starting point bias:
    {'subj_idx':3, 'v':v, 'dc':dc, 'z':0.52*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':v, 'dc':dc, 'z':0.62*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':v, 'dc':dc, 'z':0.72*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},

    # 2 DDM drift bias:
    {'subj_idx':6, 'v':v, 'dc':dc+0.04, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':v, 'dc':dc+0.05, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':v, 'dc':dc+0.06, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},

    # 3 DDM collapsing / expanding bounds:
    {'subj_idx':9, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':-0.04, 'b1_collapse':-0.04, 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':-0.05, 'b1_collapse':-0.05, 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':-0.06, 'b1_collapse':-0.06, 'nr_trials':nr_trials},

    # 4 DDM one collapsing bound:
    {'subj_idx':12, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':-0.04, 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':-0.05, 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'b0_collapse':0, 'b1_collapse':-0.06, 'nr_trials':nr_trials},

    # 5 DDM increasing drift bias:
    {'subj_idx':15, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope + 0.15, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope + 0.20, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope + 0.25, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},

    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 3
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)

groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14], [15,16,17],]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
for i, group in enumerate(groups):
    
    # neutral:
    df = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(0)))
    mean_correct = df.correct.mean()
    mean_response = df.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    # plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]

    fig = conditional_response_plot(df, quantiles, mean_response)
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    fig = summary_plot(df, quantiles, mean_correct, mean_response)
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))
