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

from sim_tools import get_OU_traces, apply_bounds_diff_trace, _bounds
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
        b1, b0 = _bounds(a=params['a'], lower_is_0=False, tmax=tmax, dt=dt)

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

# leaky regime:
v = 1
a = 0.42 # means boundary separation of 0.84 as upper boundary becomes 0.42, and lower -0.42
λ = 4
dc = 0
z = 0
ndt = 0.1

# self-excitation regime:
a2 = .6
λ2 = -2.5


sArray = [
        
    # OU neutral:
    {'subj_idx':0, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # OU starting point bias:
    {'subj_idx':3, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.04,z], 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.12,z], 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc,dc], 'z':[z+0.20,z], 'nr_trials':nr_trials},

    # OU input bias:
    {'subj_idx':6, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.2,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.5,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'λ':[λ,λ], 'a':a, 'dc':[dc+0.8,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # OU leak bias:
    {'subj_idx':9, 'v':[v,0], 'λ':[λ-2,λ+2], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'λ':[λ-3,λ+3], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'λ':[λ-4,λ+4], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    
    # OU neutral:
    {'subj_idx':12, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # OU starting point bias:
    {'subj_idx':15, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z+0.05,z], 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z+0.10,z], 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc,dc], 'z':[z+0.15,z], 'nr_trials':nr_trials},

    # OU input bias:
    {'subj_idx':18, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc+0.2,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc+0.5,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':[v,0], 'λ':[λ2,λ2], 'a':a2, 'dc':[dc+0.8,dc], 'z':[z,z], 'nr_trials':nr_trials},

    # # OU leak bias:
    {'subj_idx':21, 'v':[v,0], 'λ':[λ2-0.5,λ2+0.5], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':22, 'v':[v,0], 'λ':[λ2-1.5,λ2+1.5], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':23, 'v':[v,0], 'λ':[λ2-2.5,λ2+2.5], 'a':a2, 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},

    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 6
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    # do_simulations(sArray[0])


# groups = [[0,1,2],]
groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14], [15,16,17], [18,19,20], [21,22,23],]
# groups = []
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

    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))
    
# save combined for DDM fitting:
groups = [[4], [7], [10], [16], [19], [22],]
for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.loc[:,'condition'] = 0
    df.loc[:,'subj_idx'] = 0
    df.to_csv(os.path.join(os.path.expanduser('~/Desktop/simulations/ddm_fits_data/'), '2018_ou_data_{}.csv'.format(i+1)))