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

from sim_tools import OU_traces_get, one_accumulater_traces_apply_bounds
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

# data_folder = '/home/degee/research/model_simulations/ou_data/'
# fig_folder = '/home/degee/research/model_simulations/ou_figs/'
data_folder = '/Users/janwillem/Desktop/simulations/ou_data/'
fig_folder = '/Users/janwillem/Desktop/simulations/ou_figs/'

t = 1000
dt = 1
timesteps = int(t/dt)
print(timesteps)
    
def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    for stim in [1,0]:
        
        x1 = OU_traces_get(v=params['v'],
                                ou=params['ou'],
                                dc=params['dc'],
                                z=params['z'],
                                pre_generated=False,
                                stim=stim,
                                nr_trials=params['nr_trials'],)

        rt_dum, response_dum = one_accumulater_traces_apply_bounds(x1, a=params['a'], b0_is_0=False)
        
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

simulate = True
nr_trials = int(1e4) #100K
ndt = 15

# # for noise_sd == 0.1:
# v = 0.01
# ou = -0.05
# a = 0.65
# dc = 0
# z = 0

# for noise_sd == 0.1 / 2:
v = 0.005
ou = -0.025
a = 0.325
dc = 0
z = 0

sArray = [
        
    # OU neutral:
    {'subj_idx':0, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # OU starting point bias:
    {'subj_idx':3, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z+0.025,z], 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z+0.075,z], 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc,dc], 'z':[z+0.150,z], 'nr_trials':nr_trials},

    # OU input bias:
    {'subj_idx':6, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc+0.0015,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc+0.0025,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'ou':[ou,ou], 'a':a, 'dc':[dc+0.0035,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # OU leak bias:
    {'subj_idx':9, 'v':[v,0], 'ou':[ou+0.005,ou-0.005], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'ou':[ou+0.015,ou-0.015], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'ou':[ou+0.025,ou-0.025], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 3
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)

# groups = [[0,1,2],]
groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11],]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
for i, group in enumerate(groups):
    
    # neutral:
    df_neutral = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(0)))
    mean_correct = df_neutral.correct.mean()
    mean_response = df_neutral.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    # plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]

    fig = conditional_response_plot(df, quantiles, mean_response)
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    fig = summary_plot(df, quantiles, mean_correct, mean_response)
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))
    
# save combined for DDM fitting:
groups = [[4], [7], [10],]
for i, group in enumerate(groups): 
    df_neutral = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(1)))
    df_neutral.loc[:,'condition'] = 0
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.loc[:,'condition'] = 1
    df_combined = pd.concat((df_neutral, df))
    df_combined.loc[:,'subj_idx'] = 0
    df_combined.loc[:,'rt'] = df_combined.loc[:,'rt'] / 200
    df_combined.to_csv(os.path.join(data_folder, '2018_ou_data_{}.csv'.format(i+1)))

    

