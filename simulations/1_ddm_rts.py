#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np
import scipy as sp
import matplotlib as mpl
mpl.use("TkAgg")
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd
from IPython import embed as shell

from sim_tools import get_DDM_traces, apply_bounds_diff_trace, _bounds, _bounds_collapse_linear, _bounds_collapse_hyperbolic
from sim_tools import summary_plot, conditional_response_plot
from tqdm import tqdm

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


data_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ddm_data/')
fig_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ddm_figs/')
fits_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/fits/')

simulate = False
parallel = False
nr_trials = int(1e5) #100K
nr_trials = int(1e4) #10.000
tmax = 5
dt = 0.01

v = 1
a = 1
dc = 0
dc_slope = 0
ndt = 0.1
sv = 0.5

sArray = [

    # # 0 DDM neutral
    # {'subj_idx':0, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 1 DDM starting point bias:
    {'subj_idx':1, 'v':v, 'dc':dc, 'z':0.50*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':v, 'dc':dc, 'z':0.52*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':3, 'v':v, 'dc':dc, 'z':0.54*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':v, 'dc':dc, 'z':0.56*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':v, 'dc':dc, 'z':0.58*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':6, 'v':v, 'dc':dc, 'z':0.60*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':v, 'dc':dc, 'z':0.62*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':v, 'dc':dc, 'z':0.64*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':9, 'v':v, 'dc':dc, 'z':0.66*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':v, 'dc':dc, 'z':0.68*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':v, 'dc':dc, 'z':0.70*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 2 DDM drift bias:
    {'subj_idx':12, 'v':v, 'dc':dc+0.00, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':v, 'dc':dc+0.08, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':v, 'dc':dc+0.16, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':15, 'v':v, 'dc':dc+0.24, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':v, 'dc':dc+0.32, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':v, 'dc':dc+0.40, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':18, 'v':v, 'dc':dc+0.48, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':v, 'dc':dc+0.56, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':v, 'dc':dc+0.64, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':21, 'v':v, 'dc':dc+0.72, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':22, 'v':v, 'dc':dc+0.80, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 3 DDM increasing drift bias:
    {'subj_idx':23, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+0.0, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':24, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+0.4, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':25, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+0.8, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':26, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+1.2, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':27, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+1.6, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':28, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+2.0, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':29, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+2.4, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':30, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+2.8, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':31, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+3.2, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':32, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+3.6, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':33, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope+4.0, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # 4 DDM drift criterion + hyperbolically collapsing bounds:
    {'subj_idx':34, 'v':v, 'dc':dc+0.00, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':35, 'v':v, 'dc':dc+0.08, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':36, 'v':v, 'dc':dc+0.16, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':37, 'v':v, 'dc':dc+0.24, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':38, 'v':v, 'dc':dc+0.32, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':39, 'v':v, 'dc':dc+0.40, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':40, 'v':v, 'dc':dc+0.48, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':41, 'v':v, 'dc':dc+0.56, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':42, 'v':v, 'dc':dc+0.64, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':43, 'v':v, 'dc':dc+0.72, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},
    {'subj_idx':44, 'v':v, 'dc':dc+0.80, 'z':0.5*a*2, 'a':a*2, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_hyperbolic', 'c':1, 'nr_trials':nr_trials},

    # 5 DDM drift bias, opposing starting point bias:
    {'subj_idx':45, 'v':v, 'dc':dc+0.00, 'z':0.500*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':46, 'v':v, 'dc':dc+0.08, 'z':0.495*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':47, 'v':v, 'dc':dc+0.16, 'z':0.490*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':48, 'v':v, 'dc':dc+0.24, 'z':0.485*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':49, 'v':v, 'dc':dc+0.32, 'z':0.480*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':50, 'v':v, 'dc':dc+0.40, 'z':0.475*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':51, 'v':v, 'dc':dc+0.48, 'z':0.470*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':52, 'v':v, 'dc':dc+0.56, 'z':0.465*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':53, 'v':v, 'dc':dc+0.64, 'z':0.460*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':54, 'v':v, 'dc':dc+0.72, 'z':0.455*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':55, 'v':v, 'dc':dc+0.80, 'z':0.450*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # # 4 DDM collapsing / expanding bounds:
    # {'subj_idx':12, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.2, 'c0':-0.2, 'nr_trials':nr_trials},
    # {'subj_idx':13, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.5, 'c0':-0.5, 'nr_trials':nr_trials},
    # {'subj_idx':14, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0.8, 'c0':-0.8, 'nr_trials':nr_trials},

    # # 5 DDM one collapsing bound:
    # {'subj_idx':15, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.3, 'nr_trials':nr_trials},
    # {'subj_idx':16, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.5, 'nr_trials':nr_trials},
    # {'subj_idx':17, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'dc_slope':dc_slope, 'sv':sv, 'bound':'collapse_linear', 'c1':0, 'c0':-0.7, 'nr_trials':nr_trials},

    ]

if simulate:
	if not parallel:
		for i, s in tqdm(enumerate(sArray)):
			do_simulations(s) 
	else:
	    from joblib import Parallel, delayed
	    n_jobs = 42
	    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
	    # do_simulations(sArray[0])

groups = [list(np.arange(1,12)), list(np.arange(12,23)), list(np.arange(23,34)), 
	list(np.arange(34,45)), list(np.arange(45,56))]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

cmaps = ["Greens", 'Blues', 'Oranges', 'Purples', 'RdPu']

for i, group in enumerate(groups):
    
    # neutral:
    df = pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(0)))
    mean_correct = df.correct.mean()
    mean_response = df.response.mean()
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    # plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9,]

    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7), cmap=cmaps[i])
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    # fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    # fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))

# save combined for DDM fitting:
groups = [list(np.arange(1,12)), list(np.arange(12,23)), list(np.arange(23,34)), list(np.arange(34,45)), list(np.arange(45,56))]
for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.to_csv(os.path.join(fits_folder, '2018_ddm_data_{}.csv'.format(i+1)))

# tmax = 1
# dt = 0.01

# # linear collapse:
# fig = plt.figure()
# t = np.arange(0,tmax,dt)
# b1, b0 = _bounds_collapse_linear(a=2, c1=0.3, c0=0.3, tmax=tmax, dt=dt)
# plt.plot(t,b1)
# plt.plot(t,b0)
# fig.savefig(os.path.join(fig_folder, 'collapse_linear.pdf'))
#
# # hyperbolic collapse:
# fig = plt.figure()
# t = np.arange(0,tmax,dt)
# b1, b0 = _bounds_collapse_hyperbolic(a=2, c=1, lower_is_0=True, tmax=tmax, dt=dt)
# plt.plot(t,b1)
# plt.plot(t,b0)
# fig.savefig(os.path.join(fig_folder, 'collapse_hyperbolic.pdf'))

# load ddm results:
for i, group in enumerate(groups):

    print(i)
    print(group)
    
    # simulated data:
    df = pd.read_csv(os.path.join(fits_folder, '2018_ddm_data_{}.csv'.format(i+1)))
    
    # model params:
    params = []
    for v in range(4):
        param = pd.read_csv(os.path.join(fits_folder, '2018_ddm_data_{}_{}_params_flat.csv'.format(i+1, v)))
        param['version'] = v
        params.append(param)
    param = pd.concat(params)

    param['z'] = param['z'] - 0.5
    for v in [1,2,3]:

        param.loc[param['version']==v, 'Dbic'] = np.array(param.loc[param['version']==v, 'bic']) - np.array(
            param.loc[param['version']==0, 'bic'])
        param.loc[param['version']==v, 'Daic'] = np.array(param.loc[param['version']==v, 'aic']) - np.array(
            param.loc[param['version']==0, 'aic'])
        param.loc[param['version']==v, 'Dbic_info'] = np.array(param.loc[param['version']==v, 'bic_info']) - np.array(
            param.loc[param['version']==0, 'bic_info'])

    param.groupby(['version'])['Dbic', 'Daic', 'Dbic_info'].mean()

    print('plotting DDM results')

    # plots:
    # 1. PARAMETER ESTIMATES
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    plt.axhline(0, xmin=-0.1, xmax=1.1, lw=0.5, color='k')
    sns.barplot(data=param.loc[param['version']==3,:].loc[:,['z', 'dc']], 
    	palette=['forestgreen', 'royalblue'], ci=None, ax=ax)
    for s in range(11):
        ax.scatter([0,1], np.array(param.loc[param['version']==3,:].loc[:,['z', 'dc']].iloc[s]), 
        	color=sns.color_palette("Greys",11)[s], linewidth=0.5, edgecolor=None, zorder=10)
    plt.xticks([0,1], ['z', '$\mathregular{v_{bias}}$'], fontsize='medium')
    plt.ylim([-0.15, 1.15])
    plt.ylabel('Parameter estimate (a.u.)')
    sns.despine(offset=5, trim=True)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_folder, 'bars_{}.pdf'.format(i+1)))

    palette = ['forestgreen', 'royalblue', 'darkcyan']
    # BIC COMPUTED DIFFERENTLY
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    plt.axhline(0, xmin=-0.1, xmax=1.1, lw=0.5, color='k')
    #sns.stripplot(x='version', y='bic', data=param.loc[param['version']!=0,:], color='lightgrey', linewidth=0.5, edgecolor='black', ax=ax)
    sns.barplot(x=np.arange(3), y=np.array(param.loc[param['version']!=0,:].groupby('version').mean()['Dbic_info']),
    	palette=palette, ci=None, ax=ax)

    # find the lowest BIC
    avgbic = np.array(param.loc[param['version']!=0,:].groupby('version').mean()['Dbic_info'])
    ax.bar(np.argmin(avgbic), np.min(avgbic), facecolor=palette[np.argmin(avgbic)], edgecolor="k")

    # plt.bar(np.arange(3), np.array(param.loc[param['version']!=0,:].groupby('version').mean()['bic']), where='mid', lw=1, color='k')
    plt.ylabel('$\mathregular{\Delta BIC}}$')
    plt.xticks(np.arange(3), ['z', '$\mathregular{v_{bias}}$', 'both'], fontsize='medium')
    plt.xlabel('')
    sns.despine(offset=5, trim=True)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_folder, 'bic_info_{}.pdf'.format(i+1)))

    # 3. correlations
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    if sp.stats.pearsonr(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'])[1] < 0.05:
        sns.regplot(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'], 
        	fit_reg=True, color='forestgreen', ax=ax, ci=None)
    else:
        sns.regplot(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'], 
        	fit_reg=False, color='forestgreen', ax=ax, ci=None)
    plt.ylim([0.45, 0.75])
    plt.ylabel('P(bias)')

    ax.spines['bottom'].set_color('forestgreen')
    ax.spines['top'].set_color('forestgreen')
    ax.xaxis.label.set_color('forestgreen')
    ax.tick_params(axis='x', colors='forestgreen')

    ax = ax.twiny()
    if sp.stats.pearsonr(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'])[1] < 0.05:
        sns.regplot(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'], 
        	fit_reg=True, ci=None, color='royalblue', ax=ax)
    else:
        sns.regplot(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'], 
        	fit_reg=False, ci=None, color='royalblue', ax=ax)
    plt.xlabel('$\mathregular{v_{bias}}$')

    ax.spines['bottom'].set_color('royalblue')
    ax.spines['top'].set_color('royalblue')
    ax.xaxis.label.set_color('royalblue')
    ax.tick_params(axis='x', colors='royalblue')

    # ax.ylim([0.45, 0.75])
    sns.despine(offset=2, trim=True, top=False)
    plt.tight_layout()
    fig.savefig(os.path.join(fig_folder, 'regs_{}.pdf'.format(i+1)))

