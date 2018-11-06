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

from sim_tools import get_DDM_traces, apply_bounds_diff_trace, _bounds, _bounds_collapse_linear, _bounds_collapse_hyperbolic
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
        x = get_DDM_traces(v=params['v'],
                            z=params['z'],
                            dc=params['dc'],
                            dc_slope=params['dc_slope'],
                            sv=params['sv'],
                            stim=stim,
                            nr_trials=params['nr_trials'],
                            tmax=tmax,
                            dt=dt,)
        
        # get bounds:
        if params['bound'] == 'default':
            b1, b0 = _bounds(a=params['a'], tmax=tmax, dt=dt)
        elif params['bound'] == 'collapse_linear':
            b1, b0 = _bounds_collapse_linear(a=params['a'], c1=params['c1'], c0=params['c0'], tmax=tmax, dt=dt)
        elif params['bound'] == 'collapse_hyperbolic':
            b1, b0 = _bounds_collapse_hyperbolic(a=params['a'], c=params['c'], tmax=tmax, dt=dt)
        
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

data_folder = os.path.expanduser('~/Desktop/simulations/ddm_data/')
fig_folder = os.path.expanduser('~/Desktop/simulations/ddm_figs/')

simulate = True
nr_trials = int(1e5) #100K
tmax = 5
dt = 0.01

v = 1
a = 1
dc = 0
dc_slope = 0
ndt = 0.1
sv = 0.5

sArray = [

    # 0 DDM neutral
    {'subj_idx':0, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':1, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 1 DDM starting point bias:
    {'subj_idx':3, 'v':v, 'dc':dc, 'z':0.56*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':v, 'dc':dc, 'z':0.62*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':v, 'dc':dc, 'z':0.68*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 2 DDM drift bias:
    {'subj_idx':6, 'v':v, 'dc':dc+0.2, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':v, 'dc':dc+0.5, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':v, 'dc':dc+0.8, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 3 DDM increasing drift bias:
    {'subj_idx':9, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+1.5, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+2.5, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+3.5, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 4 DDM collapsing / expanding bounds:
    {'subj_idx':12, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.2, 'c0':-0.2, 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.5, 'c0':-0.5, 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.8, 'c0':-0.8, 'nr_trials':nr_trials},

    # 5 DDM one collapsing bound:
    {'subj_idx':15, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.3, 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.5, 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.7, 'nr_trials':nr_trials},

    # 6 DDM drift criterion + hyperbolically collapsing bounds:
    {'subj_idx':18, 'v':v, 'dc':dc+0.2, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':v, 'dc':dc+0.5, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':v, 'dc':dc+0.8, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},

    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 6
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    # do_simulations(sArray[0])

groups = [[0,1,2], [3,4,5], [6,7,8], [9,10,11], [12,13,14], [15,16,17], [18,19,20],]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
for i, group in enumerate(groups):
    
    # neutral:
    df = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(0)))
    mean_correct = df.correct.mean()
    mean_response = df.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    # plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9,]

    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))

# # save combined for DDM fitting:
# groups = [[1], [4], [7], [10],]
# for i, group in enumerate(groups): 
#     df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
#     df.loc[:,'condition'] = 0
#     df.loc[:,'subj_idx'] = 0
#     df.to_csv(os.path.join(os.path.expanduser('~/Desktop/simulations/ddm_fits_data/'), '2018_ddm_data_{}.csv'.format(i+1)))

tmax = 1
dt = 0.01

# linear collapse:
fig = plt.figure()
t = np.arange(0,tmax,dt)
b1, b0 = _bounds_collapse_linear(a=2, c1=0.3, c0=0.3, tmax=tmax, dt=dt)
plt.plot(t,b1)
plt.plot(t,b0)
fig.savefig(os.path.join(fig_folder, 'collapse_linear.pdf'))

# hyperbolic collapse:
fig = plt.figure()
t = np.arange(0,tmax,dt)
b1, b0 = _bounds_collapse_hyperbolic(a=2, c=1, lower_is_0=True, tmax=tmax, dt=dt)
plt.plot(t,b1)
plt.plot(t,b0)
fig.savefig(os.path.join(fig_folder, 'collapse_hyperbolic.pdf'))