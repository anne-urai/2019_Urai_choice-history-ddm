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

from sim_tools import LCA_traces_get, two_accumulater_traces_apply_bounds
from sim_tools import summary_plot, conditional_response_plot, traces_plot

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

# data_folder = '/home/degee/research/model_simulations/lca_data/'
# fig_folder = '/home/degee/research/model_simulations/lca_figs/'
data_folder = '/Users/janwillem/Desktop/simulations/lca_data/'
fig_folder = '/Users/janwillem/Desktop/simulations/lca_figs/'

t = 1000
dt = 1
timesteps = int(t/dt)
print(timesteps)
    
def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    x1 = []
    x2 = []
    for stim in [1,0]:
        
        # compute:
        x1_dum, x2_dum = LCA_traces_get(v=params['v'],
                                k=params['k'],
                                w=params['w'],
                                dc=params['dc'],
                                z=params['z'],
                                pre_generated=False,
                                stim=stim,
                                nr_trials=params['nr_trials'],)
        rt_dum, response_dum = two_accumulater_traces_apply_bounds(x1_dum, x2_dum, a=params['a'],)

        # append:
        x1.append(x1_dum)
        x2.append(x2_dum)
        rt.append(rt_dum)
        response.append(response_dum)
        stimulus.append(np.ones(params['nr_trials']) * stim)
    
    # traces:
    x1 = pd.DataFrame(np.array(np.vstack(x1), dtype='float16')).iloc[:200,:]
    x2 = pd.DataFrame(np.array(np.vstack(x2), dtype='float16')).iloc[:200,:]
    if os.path.exists(os.path.join(data_folder, 'x1_{}.hdf'.format(params['subj_idx']))):
        os.remove(os.path.join(data_folder, 'x1_{}.hdf'.format(params['subj_idx'])))
    if os.path.exists(os.path.join(data_folder, 'x2_{}.hdf'.format(params['subj_idx']))):
        os.remove(os.path.join(data_folder, 'x2_{}.hdf'.format(params['subj_idx'])))
    x1.to_hdf(os.path.join(data_folder, 'x1_{}.hdf'.format(params['subj_idx'])), key='lca')
    x2.to_hdf(os.path.join(data_folder, 'x2_{}.hdf'.format(params['subj_idx'])), key='lca')

    # rt dataframe:    
    df = pd.DataFrame()
    df.loc[:,'rt'] = (np.concatenate(rt) * dt) + ndt
    df.loc[:,'response'] = np.concatenate(response)
    df.loc[:,'stimulus'] = np.concatenate(stimulus)
    df.loc[:,'correct'] = np.array(np.concatenate(stimulus) == np.concatenate(response), dtype=int)
    df.loc[:,'subj_idx'] = params['subj_idx']
    df.to_csv(os.path.join(data_folder, 'df_{}.csv'.format(params['subj_idx'])))

simulate = True
nr_trials = int(1e5) #100K
ndt = 15

v = 0.01
k = 0.05
w = 0.054 # neutral
# w = 0.054 + 0.040 # inhibition dominant
a = 1.5
dc = 0.1
z = 0

sArray = [
    
    # LCA neutral:
    {'subj_idx':0, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # LCA starting point bias:
    {'subj_idx':3, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z+0.1,z], 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z+0.2,z], 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z+0.3,z], 'nr_trials':nr_trials},

    # LCA input bias:
    {'subj_idx':6, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.004,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.005,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.006,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # LCA leak bias:
    {'subj_idx':9, 'v':[v,0], 'k':[k-0.002,k+0.002], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'k':[k-0.003,k+0.003], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'k':[k-0.004,k+0.004], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # LCA inhibition bias:
    {'subj_idx':12, 'v':[v,0], 'k':[k,k], 'w':[w+0.002,w-0.002], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':[v,0], 'k':[k,k], 'w':[w+0.003,w-0.003], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':[v,0], 'k':[k,k], 'w':[w+0.004,w-0.004], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 3
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    # for params in sArray:
    #     do_simulations(params)

groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14],]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
for i, group in enumerate(groups):
    
    # neutral:
    df_neutral = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(0)))
    mean_correct = df_neutral.correct.mean()
    mean_response = df_neutral.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    x1 = pd.concat([pd.read_hdf(os.path.join(data_folder, 'x1_{}.hdf'.format(g))) for g in group], axis=0)
    x2 = pd.concat([pd.read_hdf(os.path.join(data_folder, 'x2_{}.hdf'.format(g))) for g in group], axis=0)

    # plot summary:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    fig = conditional_response_plot(df, quantiles, mean_response)
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    # plot conditional response function:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    fig = summary_plot(df, quantiles, mean_correct, mean_response)
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))

    # plot traces:
    fig = traces_plot(df, x1, x2, a, ndt)
    fig.savefig(os.path.join(fig_folder, 'traces_{}.pdf'.format(i)))

# save combined for DDM fitting:
groups = [[4], [7], [10], [13],]
for i, group in enumerate(groups): 
    df_neutral = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(1)))
    df_neutral.loc[:,'condition'] = 0
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.loc[:,'condition'] = 1
    df_combined = pd.concat((df_neutral, df))
    df_combined.loc[:,'subj_idx'] = 0
    df_combined.loc[:,'rt'] = df_combined.loc[:,'rt'] / 200
    df_combined.to_csv(os.path.join(data_folder, '2018_lca_data_{}.csv'.format(i+1)))
