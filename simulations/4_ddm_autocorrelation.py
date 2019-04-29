#!/usr/bin/env python
# encoding: utf-8

import os
import numpy as np
import scipy as sp

import matplotlib as mpl
mpl.use("TkAgg")
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib import pyplot as plt

import seaborn as sns
import pandas as pd
from IPython import embed as shell

from datetime import datetime
from tqdm import tqdm
import random

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


# https://stackoverflow.com/questions/33898665/python-generate-array-of-specific-autocorrelation
# GENERATE A PROCESS WITH A SPECIFIC AUTOCORRELATION
def sample_signal(n_samples, corr, mu=0, sigma=1):

    # McKinney, Perktold, Seabold (statsmodels)Python Time Series AnalysisSciPy Conference 201114 / 29
    import statsmodels.tsa.arima_process as ap
    import statsmodels.tsa.stattools as stattools

    # generate a process with specific autocorrelation at lag 1
    ma_coef = [1, -.5]
    y = ap.arma_generate_sample([1, -corr], [1,0], n_samples, sigma=sigma) + mu

    # test if this worked correctly
    acf = stattools.acf(y)
    assert(np.isclose(acf[1], corr, rtol=1e-01, atol=1e-01, equal_nan=False))

    return y

def do_simulations(params):

    # SIMULATE A SLOWLY AUTOCORRELATED DRIFT BIAS
    v_autocorr  = sample_signal(params['nr_trials'], params['v_autocorr'], mu=params['v'], sigma=1)

    # SIMULATE A SEQUENCE OF RANDOM STIMULI
    stims       = np.repeat([0, 1], int(params['nr_trials']/2))
    np.random.shuffle(stims)

    ## NOW DRAW RESPONSES FROM THAT PROCESS!
    rt = []
    response = []
    stimulus = []
    for t in tqdm(range(params['nr_trials'])):
        
        # get traces:
        x = get_DDM_traces(v=v_autocorr[t],
                            z=params['z'],
                            dc=params['dc'],
                            dc_slope=0,
                            sv=params['sv'],
                            stim=stims[t],
                            nr_trials=1,
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
        
        # store results for single rtial
        rt.append((rt_dum*dt)+ndt)
        response.append(response_dum)
        stimulus.append(stims[t])

    df = pd.DataFrame()
    df.loc[:,'rt'] = np.concatenate(rt)
    df.loc[:,'response'] = np.concatenate(response)
    df.loc[:,'stimulus'] = stimulus
    df.loc[:,'correct'] = np.array(stimulus == np.concatenate(response), dtype=int)
    df.loc[:,'subj_idx'] = params['subj_idx']

    # NOW RECODE REPEATS/ALTERNATES INTO RESPONSES
    df['prevresp'] = np.roll(df['response'], 1)
    df['repeat']   = 1 * (df['prevresp'] == df['response'])
    df['choice']   = df['response'].copy()
    df['response'] = df['repeat']

    df.to_csv(os.path.join(data_folder, 'df_{}.csv'.format(params['subj_idx'])))

data_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ddm_autocorr_data/')
fig_folder  = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/ddm_autocorr_figs/')
fits_folder = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/fits/')

simulate = True
parallel = False
nr_trials = int(1e3) #100K
tmax = 5
dt = 0.01

v = 1
a = 1
dc = 0
dc_slope = 0
ndt = 0.1
sv = 0.5

sArray = [

    {'subj_idx':0, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # positive autocorrelation
    {'subj_idx':1, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.1, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':2, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.2, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':3, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.3, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':4, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.4, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':5, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.5, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':6, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.6, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':7, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.7, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':8, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.8, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':9, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':0.9, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':10, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':1, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

    # negative autocorrelation, i.e. 'regression to the mean'
    {'subj_idx':11, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.1, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':12, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.2, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.3, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.4, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':15, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.5, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':16, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.6, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':17, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.7, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':18, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.8, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-0.9, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':v, 'dc':dc, 'z':0.5*a, 'a':a, 'v_autocorr':-1, 'sv':sv, 'bound':'default', 'nr_trials':nr_trials},

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

groups = [list(np.arange(0,1)), list(np.arange(1,11)), list(np.arange(11,21))]
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
# groups = [list(np.arange(1,12)), list(np.arange(12,23)), list(np.arange(23,34)), list(np.arange(34,45)), list(np.arange(45,56))]

for i, group in enumerate(groups): 
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    df.to_csv(os.path.join(fits_folder, '2018_ddm_autocorr_data_{}.csv'.format(i+1)))

## DDMS ARE FIT HERE??




tmax = 1
dt = 0.01

# load ddm results:
for i, group in enumerate(groups):
    
    # simulated data:
    df = pd.read_csv(os.path.join(fits_folder, '2018_ddm_autocorr_data_{}.csv'.format(i+1)))
    
    # model params:
    params = []
    for v in range(4):
        param = pd.read_csv(os.path.join(fits_folder, '2018_ddm_autocorr_data_{}_{}'.format(i+1, v), 'results.csv'))
        param['version'] = v
        params.append(param)
    param = pd.concat(params)
    param['z'] = param['z'] - 0.5
    for v in [1,2,3]:
        param.loc[param['version']==v, 'bic'] = np.array(param.loc[param['version']==v, 'bic']) - np.array(param.loc[param['version']==0, 'bic'])

    # plots:
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    sns.barplot(data=param.loc[param['version']==3,:].loc[:,['z', 'dc']], palette=['forestgreen', 'royalblue'], ci=None, ax=ax)
    for s in range(11):
        ax.scatter([0,1], np.array(param.loc[param['version']==3,:].loc[:,['z', 'dc']].iloc[s]), color=sns.color_palette("Greys",11)[s], linewidth=0.5, edgecolor='black', zorder=10)
    plt.xticks([0,1], ['z_bias', 'v_bias'])
    plt.tight_layout()
    sns.despine(offset=2, trim=False)
    fig.savefig(os.path.join(fig_folder, 'bars_{}.pdf'.format(i+1)))

    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    # sns.stripplot(x='version', y='bic', data=param.loc[param['version']!=0,:], color='lightgrey', linewidth=0.5, edgecolor='black', ax=ax)
    plt.step(np.arange(3), np.array(param.loc[param['version']!=0,:].groupby('version').mean()['bic']), where='mid', lw=1, color='k')
    plt.ylabel('delta AIC')
    plt.tight_layout()
    sns.despine(offset=2, trim=False)
    fig.savefig(os.path.join(fig_folder, 'bics_{}.pdf'.format(i+1)))

    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(111)
    if sp.stats.pearsonr(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'])[1] < 0.05:
        sns.regplot(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'], fit_reg=True, color='forestgreen', ax=ax)
    else:
        sns.regplot(x=param.loc[param['version']==3,'z'], y=df.groupby('subj_idx').mean()['response'], fit_reg=False, color='forestgreen', ax=ax)
    plt.ylabel('P(bias)')
    ax = ax.twiny()
    if sp.stats.pearsonr(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'])[1] < 0.05:
        sns.regplot(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'], fit_reg=True, color='royalblue', ax=ax)
    else:
        sns.regplot(x=param.loc[param['version']==3,'dc'], y=df.groupby('subj_idx').mean()['response'], fit_reg=False, color='royalblue', ax=ax)
    plt.tight_layout()
    sns.despine(offset=2, trim=False, top=False)
    fig.savefig(os.path.join(fig_folder, 'regs_{}.pdf'.format(i+1)))

