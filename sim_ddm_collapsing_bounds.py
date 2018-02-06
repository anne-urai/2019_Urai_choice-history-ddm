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
import bottleneck as bn
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

def DDM(v=1, a=1, z=0.5, dc=0, stim=0, nr_trials=1000, dt=0.01, verbose=False):
    
    """
    DDM
    """
    
    # Setup all variables:
    rt = np.zeros(nr_trials)
    response = np.zeros(nr_trials)
    
    # Run the traces:
    for t in range(nr_trials):
        time = 0
        dv = z * a
        # trace = []
        while True:
            if verbose:
                print('dt: {}'.format(dt))
            
            # update dv:
            noise = np.random.normal(0,1) / np.sqrt(0.01)
            if stim == 1:
                delta = (v + dc + noise) * dt
            else:
                delta = (-v + dc + noise) * dt
            dv += delta
            
            # Check if one of the thresholds is crossed:
            if dv >= a:
                rt[t] = time
                response[t] = 1
                break
            elif dv <= 0:
                rt[t] = time
                response[t] = 0
                break
            
            # update time:
            time += 1
            
    return(rt, response)
def DDM2(v=1, a=1, z=0.5, dc=0, stim=0, nr_trials=1000, dt=0.01, verbose=False):
    
    """
    DDM
    """
    
    # Setup all variables:
    rt = np.zeros(nr_trials)
    response = np.zeros(nr_trials)
    
    # Run the traces:
    for t in range(nr_trials):
        time = 0
        dv = z * a
        bound1 = a
        bound2 = 0
        # trace = []
        while True:
            if verbose:
                print('dt: {}'.format(dt))
            
            # update dv:
            noise = np.random.normal(0,1) / np.sqrt(0.01)
            if stim == 1:
                delta = (v + noise) * dt
            else:
                delta = (-v + noise) * dt
            dv += delta
            
            # update bounds:
            bound1 = bound1 - (dc * dt)
            bound2 = bound2 - (dc * dt)
            
            # Check if one of the thresholds is crossed:
            if dv >= bound1:
                rt[t] = time
                response[t] = 1
                break
            elif dv <= bound2:
                rt[t] = time
                response[t] = 0
                break
            
            # update time:
            time += 1
            
    return(rt, response)
def DDM3(v=1, a=1, z=0.5, dc=0, stim=0, nr_trials=1000, dt=0.01, verbose=False):
    
    """
    DDM
    """
    
    # Setup all variables:
    rt = np.zeros(nr_trials)
    response = np.zeros(nr_trials)
    
    # Run the traces:
    for t in range(nr_trials):
        time = 0
        dv = z * a
        bound1 = a
        bound2 = 0
        # trace = []
        while True:
            if verbose:
                print('dt: {}'.format(dt))
            
            # update dv:
            noise = np.random.normal(0,1) / np.sqrt(0.01)
            if stim == 1:
                delta = (v + noise) * dt
            else:
                delta = (-v + noise) * dt
            dv += delta
            
            # update bounds:
            bound1 = bound1 - (dc * dt)
            # bound2 = bound2 - (dc * dt)
            
            # Check if one of the thresholds is crossed:
            if dv >= bound1:
                rt[t] = time
                response[t] = 1
                break
            elif dv <= bound2:
                rt[t] = time
                response[t] = 0
                break
            
            # update time:
            time += 1
            
    return(rt, response)

sArray = [
    {'v':0.7, 'dc': 0, 'z':0.5, 'a':1.8, 'nr_trials':100000},
    {'v':0.7, 'dc': 0, 'z':0.625, 'a':1.8, 'nr_trials':10000},
    {'v':0.7, 'dc': 0.3, 'z':0.5, 'a':1.8, 'nr_trials':100000},
    {'v':0.7, 'dc': 0.3, 'z':0.5, 'a':1.8, 'nr_trials':100000},
    {'v':0.7, 'dc': 0.3, 'z':0.5, 'a':1.8, 'nr_trials':100000},
    ]

dfs = []
for s in range(len(sArray)):
    print(s)
    rt = []
    response = []
    stimulus = []
    for stim in [0,1]:
        
        if s < 3:
            rt_dum, response_dum = DDM.DDM(v=sArray[s]['v'], a=sArray[s]['a'], z=sArray[s]['z'], dc=sArray[s]['dc'], stim=stim, nr_trials=sArray[s]['nr_trials'], verbose=False)
        elif s == 3:
            rt_dum, response_dum = DDM.DDM2(v=sArray[s]['v'], a=sArray[s]['a'], z=sArray[s]['z'], dc=sArray[s]['dc'], stim=stim, nr_trials=sArray[s]['nr_trials'], verbose=False)
        elif s == 4:
            rt_dum, response_dum = DDM.DDM3(v=sArray[s]['v'], a=sArray[s]['a'], z=sArray[s]['z'], dc=sArray[s]['dc'], stim=stim, nr_trials=sArray[s]['nr_trials'], verbose=False)
        rt.append(rt_dum)
        response.append(response_dum)
        stimulus.append(np.ones(sArray[s]['nr_trials']) * stim)
        
    df = pd.DataFrame()
    df.loc[:,'rt'] = (np.concatenate(rt) + 200) / 1000.0
    df.loc[:,'response'] = np.concatenate(response)
    df.loc[:,'stimulus'] = np.concatenate(stimulus)
    df.loc[:,'correct'] = np.array(np.concatenate(stimulus) == np.concatenate(response), dtype=int)
    df.loc[:,'subj_idx'] = s
    dfs.append(df)
df = pd.concat(dfs)
df.to_csv('df.csv')

fig = plt.figure(figsize=(8,6))
plt_nr = 1
for s in range(len(sArray)):
    
    data_subj = df.query('subj_idx == {}'.format(s)).copy()
    
    # rt distributions:
    ax = fig.add_subplot(3,5,plt_nr)
    ax.hist(data_subj.loc[(data_subj.response==0), 'rt'], alpha=0.5, bins=40)
    ax.hist(data_subj.loc[(data_subj.response==1), 'rt'], alpha=0.5, bins=40)
    ax.set_xlim(0,0.8)
    ax.set_title('choice={}; correct={}'.format(round(data_subj.loc[:, 'response'].mean(), 3), round(data_subj.loc[:, 'correct'].mean(), 3)))
    ax.set_xlabel('RT (s)')
    plt_nr += 1
    
for s in range(len(sArray)):
    
    data_subj = df.query('subj_idx == {}'.format(s)).copy()
    
    # condition accuracy plots:
    quantiles = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
    ax = fig.add_subplot(3,5,plt_nr)
    plt.axhline(df.loc[df.subj_idx==0, 'correct'].mean(), lw=0.5, color='k')
    data_subj.loc[:,'rt_bin'] = pd.qcut(data_subj['rt'], quantiles, labels=False)
    ax.plot(np.array(data_subj.groupby('rt_bin').mean()['rt']), np.array(data_subj.groupby('rt_bin').mean()['correct']))
    ax.set_xlim(0.2,0.45)
    ax.set_ylim(0.50, 1)
    ax.set_title('Conditional accuracy')
    ax.set_xlabel('RT (s)')
    plt_nr += 1
    
for s in range(len(sArray)):

    data_subj = df.query('subj_idx == {}'.format(s)).copy()
    
    # condition accuracy plots:
    ax = fig.add_subplot(3,5,plt_nr)
    data_subj.loc[:,'rt_bin'] = pd.qcut(data_subj['rt'], quantiles, labels=False)
    plt.axhline(df.loc[df.subj_idx==0, 'response'].mean(), lw=0.5, color='k')
    ax.plot(np.array(data_subj.groupby('rt_bin').mean()['rt']), np.array(data_subj.groupby('rt_bin').mean()['response']))
    ax.set_xlim(0.2,0.45)
    ax.set_ylim(0.25,0.75)
    ax.set_title('Conditional response')
    ax.set_xlabel('RT (s)')
    plt_nr += 1
    
sns.despine(offset=10, trim=True)
plt.tight_layout()
fig.savefig('rt_dists.pdf')