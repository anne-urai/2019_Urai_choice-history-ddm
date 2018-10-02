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

from sim_tools import get_LCA_traces, apply_bounds_accumulater_traces
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

def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    x1 = []
    x2 = []
    for stim in [1,0]:
        
        # get traces:
        x1_dum, x2_dum = get_LCA_traces(v=params['v'],
                                k=params['k'],
                                w=params['w'],
                                dc=params['dc'],
                                z=params['z'],
                                pre_generated=False,
                                stim=stim,
                                nr_trials=params['nr_trials'],
                                tmax=tmax,
                                dt=dt,)
        
        # apply bounds:
        rt_dum, response_dum = apply_bounds_accumulater_traces(x1_dum, x2_dum, a=params['a'],)

        # store results:
        x1.append(x1_dum)
        x2.append(x2_dum)
        rt.append(rt_dum)
        response.append(response_dum)
        stimulus.append(np.ones(params['nr_trials']) * stim)

    # traces:
    x1 = pd.DataFrame(np.array(np.vstack(x1), dtype='float16')).iloc[:200,:]
    x2 = pd.DataFrame(np.array(np.vstack(x2), dtype='float16')).iloc[:200,:]
    x1.columns = np.arange(0,tmax,dt)
    x2.columns = np.arange(0,tmax,dt)
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

data_folder = os.path.expanduser('~/Desktop/simulations/lca_data/')
fig_folder = os.path.expanduser('~/Desktop/simulations/lca_figs/')
   
simulate = True
nr_trials = int(1e5) #100K
tmax = 5
dt = 0.01

v = 1.4
a = 1.4
k = 5
w = 5.4 # leak dominant
dc = 10
z = 0
ndt = 0.1

a2 = 1.4 # inhibition dominant
w2 = 9.4 # inhibition dominant

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
    {'subj_idx':6, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.3,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.6,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'k':[k,k], 'w':[w,w], 'a':[a,a], 'dc':[dc+0.9,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # LCA leak bias:
    {'subj_idx':9, 'v':[v,0], 'k':[k-0.2,k+0.2], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'k':[k-0.4,k+0.4], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'k':[k-0.6,k+0.6], 'w':[w,w], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # LCA inhibition bias:
    {'subj_idx':12, 'v':[v,0], 'k':[k,k], 'w':[w+0.2,w-0.2], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':[v,0], 'k':[k,k], 'w':[w+0.4,w-0.4], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':[v,0], 'k':[k,k], 'w':[w+0.6,w-0.6], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # # LCA neutral:
    # {'subj_idx':15, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':16, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':17, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # # LCA starting point bias:
    # {'subj_idx':18, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z+0.05,z], 'nr_trials':nr_trials},
    # {'subj_idx':19, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z+0.15,z], 'nr_trials':nr_trials},
    # {'subj_idx':20, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z+0.25,z], 'nr_trials':nr_trials},
    
    # # LCA input bias:
    # {'subj_idx':21, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc+0.3,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':22, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc+0.6,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':23, 'v':[v,0], 'k':[k,k], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc+0.9,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # # LCA leak bias:
    # {'subj_idx':24, 'v':[v,0], 'k':[k-0.2,k+0.2], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':25, 'v':[v,0], 'k':[k-0.5,k+0.5], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':26, 'v':[v,0], 'k':[k-0.8,k+0.8], 'w':[w2,w2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # # LCA inhibition bias:
    # {'subj_idx':27, 'v':[v,0], 'k':[k,k], 'w':[w2+0.2,w2-0.2], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':28, 'v':[v,0], 'k':[k,k], 'w':[w2+0.5,w2-0.5], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    # {'subj_idx':29, 'v':[v,0], 'k':[k,k], 'w':[w2+0.8,w2-0.8], 'a':[a2,a2], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 6
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    # for params in sArray:
    #     do_simulations(params)

groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14],[15,16,17], [18,19,20], [21,22,23], [24,25,26], [27,28,29],]
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
    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    # plot conditional response function:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]
    fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))

    # plot traces:
    fig = traces_plot(df, x1, x2, a, ndt)
    fig.savefig(os.path.join(fig_folder, 'traces_{}.pdf'.format(i)))

# save combined for DDM fitting:
groups = [[4], [7], [10], [13], [19], [22], [25], [28],]
for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.loc[:,'condition'] = 0
    df.loc[:,'subj_idx'] = 0
    df.to_csv(os.path.join(os.path.expanduser('~/Desktop/simulations/ddm_fits_data/'), '2018_lca_data_{}.csv'.format(i+1)))
