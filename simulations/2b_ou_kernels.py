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

from models import OU_traces_get, two_accumulater_traces_apply_timepoint

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

data_folder = '/home/degee/research/model_simulations/ou_data/'
fig_folder = '/home/degee/research/model_simulations/ou_figs/'

t = 200
dt = 1
timesteps = int(t/dt)
print(timesteps)

def do_simulations(params):
    
    rt = []
    response = []
    stimulus = []
    for stim in [1,0]:
        
        x1, x2 = OU_traces_get(v=params['v'],
                                ou=params['ou'],
                                dc=params['dc'],
                                z=params['z'],
                                pre_generated=True,
                                stim=stim,
                                nr_trials=params['nr_trials'],
                                timesteps=timesteps,
                                )
        rt_dum, response_dum = two_accumulater_traces_apply_timepoint(x1, x2)
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
nr_trials = 25000
ndt = 0

v = 0
ou = -0.005
a = 0.15
dc = 0.1
z = 0

# inputs = [np.array(np.vstack([np.repeat(v, timesteps) + np.random.normal(0,0.1,timesteps) for _ in range(int(nr_trials))]), dtype='float16'),
#             np.array(np.vstack([np.repeat(0, timesteps) + np.random.normal(0,0.1,timesteps) for _ in range(int(nr_trials))]), dtype='float16'),]
# np.save(os.path.join(data_folder, 'inputs.npy'), inputs)
inputs = np.load(os.path.join(data_folder, 'inputs.npy'))

sArray = [
    # LCA:
    {'subj_idx':12, 'v':inputs, 'ou':[ou,ou], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':13, 'v':inputs, 'ou':[ou,ou], 'a':[a,a], 'dc':[dc,dc], 'z':[z+0.5,z], 'nr_trials':nr_trials},
    {'subj_idx':14, 'v':inputs, 'ou':[ou,ou], 'a':[a,a], 'dc':[dc+0.002,dc], 'z':[z,z], 'nr_trials':nr_trials},
    {'subj_idx':15, 'v':inputs, 'ou':[ou+0.0001,ou-0.0002], 'a':[a,a], 'dc':[dc,dc], 'z':[z,z], 'nr_trials':nr_trials},
    ]

if simulate:
    from joblib import Parallel, delayed
    n_jobs = 8
    res = Parallel(n_jobs=n_jobs)(delayed(do_simulations)(params) for params in sArray)
    
groups = [[12], [13], [14], [15], ]
# groups = [[0], [1], [2], ]

for i, group in enumerate(groups):
    
    subj = sArray[i]['subj_idx']
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join(data_folder, 'df_{}.csv'.format(g))) for g in group], axis=0)
    yes = np.array(df.loc[(df["subj_idx"]==subj), "response"] == 1)
    no = np.array(df.loc[(df["subj_idx"]==subj), "response"] == 0)
    stimulus = np.array(df.loc[(df["subj_idx"]==subj), "stimulus"] == 1)
    
    fluctuations = np.vstack((inputs[0]-inputs[1],inputs[1]-inputs[0]))
    kernel_yes = fluctuations[yes,:].mean(axis=0)[:-1]
    kernel_no = fluctuations[no,:].mean(axis=0)[:-1]
    
    if i == 0:
        kernel_yes_ = kernel_yes.copy()
        kernel_no_ = kernel_no.copy()
    
    x = np.arange(kernel_yes.shape[0]) * dt
    
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(1,1,1)
    if i == 0:
        ax.plot(x, kernel_yes_, lw=1, color='orange')
        ax.plot(x, kernel_no_, lw=1, color='forestgreen')
    else:
        # ax.plot(x, kernel_yes, lw=1, color='orange', label='choice a')
        # ax.plot(x, kernel_no, lw=1, color='forestgreen', label='choice b')
        ax.plot(x, kernel_yes-kernel_yes_, lw=1, color='orange', label='choice a')
        ax.plot(x, kernel_no-kernel_no_, lw=1, color='forestgreen', label='choice b')
    
    plt.axhline(0, color='k', lw=0.5)
    # ax.set_ylim(-0.25, 0.25)
    ax.set_xlim(0, t)
    ax.set_title('choice={}; correct={}'.format(round(df.loc[:, 'response'].mean(), 3), 
                                                round(df.loc[:, 'correct'].mean(), 3),))
    ax.set_xlabel('Timesteps')
    ax.set_ylabel('mean input')
    sns.despine(offset=10, trim=True)
    plt.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(fig_folder, 'kernels_{}.pdf'.format(i)))
