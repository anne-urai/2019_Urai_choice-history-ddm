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

from sim_tools import get_OU_traces, apply_bounds_diff_trace, _bounds, _bounds_collapse_linear, _bounds_collapse_hyperbolic
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
    
def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    for stim in [1,0]:
        
        # get traces:
        x = get_OU_traces(v=params['v'],
                                λ=params['λ'],
                                dc=params['dc'],
                                z=params['z'],
                                pre_generated=False,
                                stim=stim,
                                nr_trials=params['nr_trials'],
                                tmax=tmax,
                                dt=dt,)

        # get bounds:
        if params['bound'] == 'default':
            b1, b0 = _bounds(a=params['a'], lower_is_0=False, tmax=tmax, dt=dt)
        elif params['bound'] == 'collapse_linear':
            b1, b0 = _bounds_collapse_linear(a=params['a'], c1=params['c1'], c0=params['c0'], lower_is_0=False, tmax=tmax, dt=dt)
        elif params['bound'] == 'collapse_hyperbolic':
            b1, b0 = _bounds_collapse_hyperbolic(a=params['a'], c=params['c'], lower_is_0=False, tmax=tmax, dt=dt)

        # apply bounds:
        rt_dum, response_dum = apply_bounds_diff_trace(x=x, b1=b1, b0=b0)
        
        # store results:
        rt.append((rt_dum*dt)+ndt)
        response.append(response_dum)
        stimulus.append(np.ones(params['nr_trials']) * stim)
        
    df = pd.DataFrame()
    df.loc[:,'rt'] = np.concatenate(rt)
    df.loc[:,'response'] = np.concatenate(response)
    df.loc[:,'stimulus'] = np.concatenate(stimulus)
    df.loc[:,'correct'] = np.array(np.concatenate(stimulus) == np.concatenate(response), dtype=int)
    df.loc[:,'subj_idx'] = params['subj_idx']
    df.to_csv(os.path.join(data_folder, 'df_{}.csv'.format(params['subj_idx'])))

data_folder = os.path.expanduser('~/Desktop/simulations/ou_data/')
fig_folder = os.path.expanduser('~/Desktop/simulations/ou_figs/')

simulate = True
nr_trials = int(1e5) #100K
tmax = 5
dt = 0.01

# self-excitation regime:
v = 1
a = 1.2
λ = -2.5
dc = 0
z = 0
ndt = 0.1

c = 0.3

# # leaky regime:
# a = 0.42 # means boundary separation of 0.84 as upper boundary becomes 0.42, and lower -0.42
# λ = 4

sArray = [
        
    # OU neutral:
    {'subj_idx':0, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # OU starting point bias:
    {'subj_idx':3, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.05,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.10,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.15,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # OU input bias:
    {'subj_idx':6, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.2,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.5,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.8,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # # OU leak bias:
    {'subj_idx':9, 'v':[v,0], 'λ':[λ-0.5,λ+0.5], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'λ':[λ-1.5,λ+1.5], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'λ':[λ-2.5,λ+2.5], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 6
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    # do_simulations(sArray[12])

# groups = [[0,1,2],]
groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11],]
# groups = []
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
for i, group in enumerate(groups):
    
    # neutral:
    df_neutral = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(1)))

    print(df_neutral['rt'].median())

    mean_correct = df_neutral.correct.mean()
    mean_response = df_neutral.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    # plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9]

    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))
    
# save combined for DDM fitting:
groups = [[4], [7], [10],]
for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.loc[:,'condition'] = 0
    df.loc[:,'subj_idx'] = 0
    df.to_csv(os.path.join(os.path.expanduser('~/Desktop/simulations/ddm_fits_data/'), '2018_ou_data_{}.csv'.format(i+1)))

# hyperbolic collapse:
fig = plt.figure()
t = np.arange(0,tmax,dt)
b1, b0 = _bounds_collapse_hyperbolic(a=a, c=c, lower_is_0=False, tmax=tmax, dt=dt)
plt.plot(t,b1)
plt.plot(t,b0)
fig.savefig(os.path.join(fig_folder, 'collapse_hyperbolic.pdf'))