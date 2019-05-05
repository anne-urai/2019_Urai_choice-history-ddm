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
from tqdm import tqdm

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
                                ll=params['ll'],
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

data_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ou_data/')
fig_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ou_figs/')
fits_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/fits/')

simulate = False
parallel = False
nr_trials = int(1e5) #100K
nr_trials = int(1e4)
tmax = 5
dt = 0.01

# self-excitation regime:
v = 1
a = 1.2
ll = -2.5
dc = 0
z = 0
ndt = 0.1
c = 0.3

# # leaky regime:
# a = 0.42 # means boundary separation of 0.84 as upper boundary becomes 0.42, and lower -0.42
# ll = 4

sArray = [
        
    # OU neutral:
    {'subj_idx':0, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # OU starting point bias:
    {'subj_idx':1, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.00,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.02,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':3, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.04,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.06,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.08,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':6, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.10,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.12,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.14,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':9, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.16,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.18,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':11, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc,dc], 'z':[z+0.20,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # OU input bias:
    {'subj_idx':12, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.0,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.1,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.2,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':15, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.3,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.4,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.5,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':18, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.6,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.7,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.8,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':21, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+0.9,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':22, 'v':[v,0], 'll':[ll,ll], 'a':a, 'dc':[dc+1.0,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    # # OU leak bias:
    {'subj_idx':23, 'v':[v,0], 'll':[ll-0.0,ll+0.0], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':24, 'v':[v,0], 'll':[ll-0.3,ll+0.3], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':25, 'v':[v,0], 'll':[ll-0.6,ll+0.6], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':26, 'v':[v,0], 'll':[ll-0.9,ll+0.9], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':27, 'v':[v,0], 'll':[ll-1.2,ll+1.2], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':28, 'v':[v,0], 'll':[ll-1.5,ll+1.5], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':29, 'v':[v,0], 'll':[ll-1.8,ll+1.8], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':30, 'v':[v,0], 'll':[ll-2.1,ll+2.1], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':31, 'v':[v,0], 'll':[ll-2.4,ll+2.4], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':32, 'v':[v,0], 'll':[ll-2.7,ll+2.7], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},
    {'subj_idx':33, 'v':[v,0], 'll':[ll-3.0,ll+3.0], 'a':a, 'dc':[dc,dc], 'z':[z,z], 'bound':'collapse_hyperbolic', 'c':c, 'nr_trials':nr_trials},

    ]

if simulate:
    if not parallel:
        for i, s in tqdm(enumerate(sArray)):
            do_simulations(s) 
    else:
        from joblib import Parallel, delayed
        n_jobs = 42
        res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
        # do_simulations(sArray[12])

# groups = [[0,1,2],]
groups = [list(np.arange(1,12)), list(np.arange(12,23)), list(np.arange(23,34)),]
quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
cmaps = ["Greens", 'Blues', 'RdPu',]
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

    fig = conditional_response_plot(df, quantiles, mean_response, xlim=(0.1,0.7), cmap=cmaps[i])
    fig.savefig(os.path.join(fig_folder, 'crf_{}.pdf'.format(i)))

    # fig = summary_plot(df, quantiles, mean_correct, mean_response, xlim=(0.1,0.7))
    # fig.savefig(os.path.join(fig_folder, 'summary_{}.pdf'.format(i)))

# save combined for DDM fitting:
groups = [list(np.arange(1,12)), list(np.arange(12,23)), list(np.arange(23,34)),]
for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.to_csv(os.path.join(fits_folder, '2018_ou_data_{}.csv'.format(i+1)))

# hyperbolic collapse:
fig = plt.figure()
t = np.arange(0,tmax,dt)
b1, b0 = _bounds_collapse_hyperbolic(a=a, c=c, lower_is_0=False, tmax=tmax, dt=dt)
plt.plot(t,b1)
plt.plot(t,b0)
fig.savefig(os.path.join(fig_folder, 'collapse_hyperbolic.pdf'))



# load ddm results:
for i, group in enumerate(groups):
    
    # simulated data:
    df = pd.read_csv(os.path.join(fits_folder, '2018_ou_data_{}.csv'.format(i+1)))
    
    # model params:
    params = []
    for v in range(4):
        param = pd.read_csv(os.path.join(fits_folder, '2018_ou_data_{}_{}_params_flat.csv'.format(i+1, v)))
        param['version'] = v
        params.append(param)

    param = pd.concat(params)
    param['z'] = param['z'] - 0.5
    for v in [1,2,3]:
        # param.loc[param['version']==v, 'bic'] = np.array(param.loc[param['version']==v, 'bic']) - np.array(param.loc[param['version']==0, 'bic'])
        # param.loc[param['version']==v, 'aic'] = np.array(param.loc[param['version']==v, 'aic']) - np.array(param.loc[param['version']==0, 'aic'])
        param.loc[param['version']==v, 'Dbic_info'] = np.array(param.loc[param['version']==v, 'bic_info']) - np.array(
            param.loc[param['version']==0, 'bic_info'])
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

    # 2. BIC
    palette = ['forestgreen', 'royalblue', 'darkcyan']
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
    fig.savefig(os.path.join(fig_folder, 'bics_{}.pdf'.format(i+1)))

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
