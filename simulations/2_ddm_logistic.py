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

def simulate_DDM(v=1, a=1, z=0.5, dc=0, dc_slope=0, b0_collapse=0, b1_collapse=0, stim=0, nr_trials=1000, dt=0.01):
    
    """
    DRIFT DIFFUSION MODEL
    
    v: drift rate
    a: boundary separation
    z: starting point
    dc: drift bias
    
    dc_slope: increase of drift bias per timeunit
    b0_collapse: collapse lower bound
    b1_collapse: collapse upper bound
    
    stim: stimulus [0 / 1]
    
    """
    
    rt = np.zeros(nr_trials)
    response = np.zeros(nr_trials)
    for t in range(nr_trials):
        time = 0
        dv = z * a
        b0 = 0
        b1 = a
        while True:
            
            # update dc:
            dc_ = (dc_slope * dt * time) + dc
            
            # update noise:
            noise = np.random.normal(0,1) * np.sqrt(dt)
            
            # updating term:
            if stim == 1:
                delta = ((v + dc_) * dt) + noise
            else:
                delta = ((-v + dc_) * dt) + noise
            
            # update dv:
            dv += delta
            
            # update bounds:
            b0 += (b0_collapse * dt)
            b1 += (b1_collapse * dt)
            
            # Check if one of the thresholds is crossed:
            if dv >= b1:
                rt[t] = time
                response[t] = 1
                break
            elif dv <= b0:
                rt[t] = time
                response[t] = 0
                break
            
            # update time:
            time += 1
            
    return(rt, response)

def simulate_DDM2(v=1, a=1, z=0.5, dc=0, dc_slope=0, b0_collapse=0, b1_collapse=0, stim=0, nr_trials=1000, dt=0.01):
    
    """
    DRIFT DIFFUSION MODEL
    
    v: drift rate
    a: boundary separation
    z: starting point
    dc: drift bias
    
    dc_slope: increase of drift bias per timeunit
    b0_collapse: collapse lower bound
    b1_collapse: collapse upper bound
    
    stim: stimulus [0 / 1]
    
    """
    
    rt = np.zeros(nr_trials)
    response = np.zeros(nr_trials)
    for t in range(nr_trials):
        time = 0
        dv = z * a
        b0 = 0
        b1 = a
        while True:
            
            # update dc:
            dc_ = (dc_slope * dt * time) + dc
            dc_scaling = 1 / (abs(v)+1)
            
            # update noise:
            noise = np.random.normal(0,1) * np.sqrt(dt)
            
            # updaring term:
            if stim == 1:
                delta = (v * dt) + (dc_scaling * dc_ * dt) + noise
            else:
                delta = (-v * dt) + (dc_scaling * dc_ * dt) + noise
            
            # update dv:
            dv += delta
            
            # update bounds:
            b0 += (b0_collapse * dt)
            b1 += (b1_collapse * dt)
            
            # Check if one of the thresholds is crossed:
            if dv >= b1:
                rt[t] = time
                response[t] = 1
                break
            elif dv <= b0:
                rt[t] = time
                response[t] = 0
                break
            
            # update time:
            time += 1
            
    return(rt, response)

# nr_trials = 20000
nr_trials = 5000
t = 0.2
sArray = [
    
    # # DDM neutral
    # {'subj_idx':0, 'v':0.0, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':1, 'v':0.3, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':2, 'v':0.6, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':3, 'v':0.9, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':4, 'v':1.2, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':5, 'v':1.5, 'dc':0, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},

    # # DDM starting point bias:
    # {'subj_idx':6, 'v':0.0, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':7, 'v':0.3, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':8, 'v':0.6, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':9, 'v':0.9, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':10, 'v':1.2, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':11, 'v':1.5, 'dc':0, 'z':0.67, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    #
    # # DDM drift bias:
    # {'subj_idx':12, 'v':0.0, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':13, 'v':0.3, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':14, 'v':0.6, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':15, 'v':0.9, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':16, 'v':1.2, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    # {'subj_idx':17, 'v':1.5, 'dc':0.4, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    
    # DDM drift bias:
    {'subj_idx':18, 'v':0.0, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':19, 'v':0.3, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':20, 'v':0.6, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':21, 'v':0.9, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':22, 'v':1.2, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    {'subj_idx':23, 'v':1.5, 'dc':0.6, 'z':0.5, 'a':1.8, 'dc_slope':0, 'b0_collapse':0, 'b1_collapse':0, 'nr_trials':nr_trials},
    
    ]

simulate = True
for s in range(len(sArray)):
    print(s)
    
    if simulate:
    
        rt = []
        response = []
        stimulus = []
        for stim in [0,1]:

            if sArray[s]['subj_idx'] <= 17:
                rt_dum, response_dum = simulate_DDM(v=sArray[s]['v'],
                                                    a=sArray[s]['a'],
                                                    z=sArray[s]['z'],
                                                    dc=sArray[s]['dc'],
                                                    dc_slope=sArray[s]['dc_slope'],
                                                    b0_collapse=sArray[s]['b0_collapse'],
                                                    b1_collapse=sArray[s]['b1_collapse'],
                                                    stim=stim,
                                                    nr_trials=sArray[s]['nr_trials'],)
        
            else:
                rt_dum, response_dum = simulate_DDM2(v=sArray[s]['v'],
                                                    a=sArray[s]['a'],
                                                    z=sArray[s]['z'],
                                                    dc=sArray[s]['dc'],
                                                    dc_slope=sArray[s]['dc_slope'],
                                                    b0_collapse=sArray[s]['b0_collapse'],
                                                    b1_collapse=sArray[s]['b1_collapse'],
                                                    stim=stim,
                                                    nr_trials=sArray[s]['nr_trials'],)
            rt.append(rt_dum)
            response.append(response_dum)
            stimulus.append(np.ones(sArray[s]['nr_trials']) * stim)

        df = pd.DataFrame()
        df.loc[:,'rt'] = (np.concatenate(rt) / 1000.0) + t
        df.loc[:,'response'] = np.concatenate(response)
        df.loc[:,'stimulus'] = np.concatenate(stimulus)
        df.loc[:,'correct'] = np.array(np.concatenate(stimulus) == np.concatenate(response), dtype=int)
        df.loc[:,'subj_idx'] = sArray[s]['subj_idx']
        df.to_csv(os.path.join('dataframes_b', 'df_{}.csv'.format(sArray[s]['subj_idx'])))


groups = [[0,1,2,3,4,5], [6,7,8,9,10,11], [12,13,14,15,16,17]]
colors = ['grey', 'lightsteelblue', 'lightcoral']
labels = ['unbiased', 'z-bias', 'dc-bias',]
fig = plt.figure(figsize=(2,4))
for i, group in enumerate(groups):
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join('dataframes_b', 'df_{}.csv'.format(g))) for g in group], axis=0)
    print(df.correct.mean())
    print(df.response.mean())
    d = df.groupby(['subj_idx', 'stimulus']).mean()
    choice = np.concatenate((d["response"][::2][::-1], d["response"][1::2]))
    correct = np.concatenate((d["correct"][::2][::-1], d["correct"][1::2]))
    v = np.array([-1.5, -1.2, -0.9, -0.6, -0.3, 0, 0, 0.3, 0.6, 0.9, 1.2, 1.5])
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(v[6:], (correct[:6][::-1]+correct[6:])/2.0, 'o', ms=3, color=colors[i], label=labels[i])
    ax.set_xlabel('Drift rate')
    ax.set_ylabel('% correct')
    ax.set_ylim(0,1)
    ax.set_xlim(0,2)
    ax.legend()
    ax = fig.add_subplot(2,1,2)
    ax.plot(v, choice, 'o', ms=3, color=colors[i], label=labels[i])
    ax.set_xlabel('Drift rate')
    ax.set_ylabel('% choice a')
    ax.set_ylim(0,1)
    ax.set_xlim(-2,2)

sns.despine(offset=10, trim=True)
plt.tight_layout()
fig.savefig('psychometric_curve.pdf')
    
groups = [[12,13,14,15,16,17], [18,19,20,21,22,23]]
colors = ['lightcoral', 'orange']
labels = ['dc-bias', 'dc-bias M',]
fig = plt.figure(figsize=(2,4))
for i, group in enumerate(groups):
    
    # load group:
    df = pd.concat([pd.read_csv(os.path.join('dataframes_b', 'df_{}.csv'.format(g))) for g in group], axis=0)
    
    print(df.correct.mean())
    print(df.response.mean())
    
    d = df.groupby(['subj_idx', 'stimulus']).mean()
    choice = np.concatenate((d["response"][::2][::-1], d["response"][1::2]))
    correct = np.concatenate((d["correct"][::2][::-1], d["correct"][1::2]))
    v = np.array([-1.5, -1.2, -0.9, -0.6, -0.3, 0, 0, 0.3, 0.6, 0.9, 1.2, 1.5])
    
    ax = fig.add_subplot(2,1,1)
    ax.plot(v[6:], (correct[:6][::-1]+correct[6:])/2.0, 'o', ms=3, color=colors[i], label=labels[i])
    ax.set_xlabel('Drift rate')
    ax.set_ylabel('% correct')
    ax.set_ylim(0,1)
    ax.set_xlim(0,2)
    ax.legend()
    ax = fig.add_subplot(2,1,2)
    ax.plot(v, choice, 'o', ms=3, color=colors[i], label=labels[i])
    ax.set_xlabel('Drift rate')
    ax.set_ylabel('% choice a')
    ax.set_ylim(0,1)
    ax.set_xlim(-2,2)

sns.despine(offset=10, trim=True)
plt.tight_layout()
fig.savefig('psychometric_curve_2.pdf')

    
