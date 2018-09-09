import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from IPython import embed as shell

def DDM_traces_get(v=1, z=0.5, dc=0, dc_slope=0, noise_sd=0.1, stim=0, nr_trials=1000, t=5.0, dt=0.01):
    
    """
    DDM

    v:  mean drift rate
    z:  starting point
    dc: drift criterion
    """
    
    if stim == 0:
        v = -v
    
    x1 = np.zeros((nr_trials, int(t/dt)))
    x1[:,:] = np.NaN
    x1[:,0] = z
    for i in range((int(t/dt))-1):
        x1[:,i+1] = x1[:,i] + ((v + dc + (dc_slope * dt * i) ) * dt) + (np.random.normal(0, noise_sd, nr_trials)*np.sqrt(dt))
    return x1
    
def OU_traces_get(v, ou, dc, z, noise_sd=0.1, pre_generated=False, stim=0, nr_trials=1000, timesteps=200,):
    
    """
    OU
    
    v:  mean drift rate
    ou: Ornstein-Uhlenbeck process parameter (effective leak / self-excitation)
    z:  starting point
    dc: drift criterion
    """
    
    if stim == 0:
        v = v[::-1]
    
    x1 = np.zeros((nr_trials, timesteps))
    x2 = np.zeros((nr_trials, timesteps))
    x1[:,:] = np.NaN
    x2[:,:] = np.NaN
    x1[:,0] = z[0]
    x2[:,0] = z[1]
    for i in range(timesteps-1):
        if pre_generated:
            x1[:,i+1] = x1[:,i] + v[0][:,i] + dc[0] + (ou[0]*x1[:,i])
            x2[:,i+1] = x2[:,i] + v[1][:,i] + dc[1] + (ou[1]*x2[:,i])
        else:
            x1[:,i+1] = x1[:,i] + v[0] + dc[0] + (ou[0]*x1[:,i]) + np.random.normal(0, noise_sd/2.0, nr_trials)
            x2[:,i+1] = x2[:,i] + v[1] + dc[1] + (ou[1]*x2[:,i]) + np.random.normal(0, noise_sd/2.0, nr_trials)
    return x1-x2

def LCA_traces_get(v, k, w, dc, z, noise_sd=0.1, pre_generated=False, stim=0, nr_trials=1000, timesteps=200,):
    
    """
    LCA
    """
    
    if stim == 0:
        v = v[::-1]
    
    x1 = np.zeros((nr_trials, timesteps))
    x2 = np.zeros((nr_trials, timesteps))
    x1[:,:] = np.NaN
    x2[:,:] = np.NaN
    x1[:,0] = z[0]
    x2[:,0] = z[1]
    for i in range(timesteps-1):
        if pre_generated:
            x1[:,i+1] = np.clip(x1[:,i] + v[0][:,i] + dc[0] - (k[0]*x1[:,i]) - (w[1]*x2[:,i]), a_min=0, a_max=1e6)
            x2[:,i+1] = np.clip(x2[:,i] + v[1][:,i] + dc[1] - (k[1]*x2[:,i]) - (w[0]*x1[:,i]), a_min=0, a_max=1e6)
        else:
            x1[:,i+1] = np.clip(x1[:,i] + v[0] + dc[0] - (k[0]*x1[:,i]) - (w[1]*x2[:,i]) + np.random.normal(0, noise_sd, nr_trials), a_min=0, a_max=1e6)
            x2[:,i+1] = np.clip(x2[:,i] + v[1] + dc[1] - (k[1]*x2[:,i]) - (w[0]*x1[:,i]) + np.random.normal(0, noise_sd, nr_trials), a_min=0, a_max=1e6)
    return x1, x2

def one_accumulater_traces_apply_bounds(x1, a=0.15, b0_is_0=True, b0_collapse=0, b1_collapse=0, dt=0.01):
    
    rt = np.zeros(x1.shape[0])
    response = np.zeros(x1.shape[0])
    for i in range(x1.shape[0]):
        bound1 = a
        if b0_is_0:
            bound0 = 0
        else:
            bound0 = -a
        for j in range(x1.shape[1]):
            if x1[i,j] >= bound1:
                rt[i] = j
                response[i] = 1
                break
            elif x1[i,j] <= bound0:
                rt[i] = j
                response[i] = 0
                break
            bound1 += b1_collapse * dt
            bound0 += b0_collapse * dt
    return rt, response

def two_accumulater_traces_apply_bounds(x1, x2, a=[0.15, 0.15],):
    
    rt = np.zeros(x1.shape[0])
    response = np.zeros(x1.shape[0])
    for i in range(x1.shape[0]):
        for j in range(x1.shape[1]):
            if x1[i,j] >= a[0]:
                rt[i] = j
                response[i] = 1
                break
            elif x2[i,j] >= a[1]:
                rt[i] = j
                response[i] = 0
                break
    return rt, response

def two_accumulater_traces_apply_timepoint(x1, x2, timepoint=None):
    
    if timepoint:
        response = np.array(x1[:,timepoint] > x2[:,timepoint], dtype=int)
        rt = np.ones(x1.shape[0]) * timepoint
    else:
        response = np.array(x1[:,-1] > x2[:,-1], dtype=int)
        rt = np.ones(x1.shape[0]) * x1.shape[1]
    
    return rt, response

def summary_plot(df, quantiles, mean_correct, mean_response):

    fig = plt.figure(figsize=(2,6))

    # rt distributions:
    ax = fig.add_subplot(3,1,1)
    ax.hist(df.loc[(df.response==0), 'rt']*-1.0, color='forestgreen', alpha=0.5, bins=25)
    ax.hist(df.loc[(df.response==1), 'rt'], color='orange', alpha=0.5, bins=25)
    ax.set_title('P(bias)={}; P(correct)={}'.format(round(df.loc[:, 'response'].mean(), 3), round(df.loc[:, 'correct'].mean(), 3)))
    ax.set_xlabel('RT (s)')
    ax.set_ylabel('Trials (#)')

    # condition accuracy plots:
    ax = fig.add_subplot(3,1,2)
    plt.axhline(mean_correct, lw=0.5, color='k')
    df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
    ax.plot(np.array(quantiles)[1:], np.array(df.groupby('rt_bin').mean()['correct']), color='k')
    plt.xticks(np.array(quantiles)[1:], list(np.array(quantiles)[1:]))
    ax.set_ylim(0.50, 1)
    ax.set_title('Conditional accuracy')
    ax.set_xlabel('RT (quantiles)')
    ax.set_ylabel('P(correct)')

    # condition response plots:
    ax = fig.add_subplot(3,1,3)
    df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
    plt.axhline(mean_response, lw=0.5, color='k')
    ax.plot(np.array(quantiles)[1:], np.array(df.groupby('rt_bin').mean()['response']), color='k')
    plt.xticks(np.array(quantiles)[1:], list(np.array(quantiles)[1:]))
    ax.set_ylim(0.25,0.75)
    ax.set_title('Conditional response')
    ax.set_xlabel('RT (quantiles)')
    ax.set_ylabel('P(bias)')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()

    return fig

def conditional_response_plot(df, quantiles, mean_response):
    
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(1,1,1)
    plt.axhline(mean_response, lw=0.5, color='k')
    df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
    d = df.groupby(['subj_idx', 'rt_bin']).mean().reset_index()
    for s, a in zip(np.unique(d["subj_idx"]), [0.1, 0.5, 0.9]):
        ax.plot(np.array(quantiles)[1:], d.loc[d["subj_idx"]==s, "response"], color='k', alpha=a)
    plt.xticks(np.array(quantiles)[1:], list(np.array(quantiles)[1:]))
    ax.set_ylim(0.25,0.75)
    ax.set_title('P(correct) = {}\nP(bias) = {}'.format(
                                                    round(df.loc[df["subj_idx"]==np.unique(df["subj_idx"])[1], 'correct'].mean(), 2),
                                                    round(df.loc[df["subj_idx"]==np.unique(df["subj_idx"])[1], 'response'].mean(), 2),
                                                    ))
    ax.set_xlabel('RT (quantiles)')
    ax.set_ylabel('P(bias)')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    
    return fig

def traces_plot(df, x1, x2, a, ndt):
    
    fig = plt.figure(figsize=(2,2))
    ax = fig.add_subplot(1,1,1)
    t = 0
    plotted = 0
    while plotted < 20:
        if df.iloc[t]['response'] == 1:
            x1_trace = x1.iloc[t]
            x2_trace = x2.iloc[t]
            cutoff = int(df.iloc[t]['rt'] - ndt)
            x1_trace.loc[x1.columns > cutoff] = np.NaN
            x2_trace.loc[x1.columns > cutoff] = np.NaN
            ax.plot(x1_trace, lw=0.66, color='r')
            ax.plot(x2_trace, lw=0.33, color='b', alpha=0.5)
            plotted += 1
        t += 1
    plt.axhline(a, lw=1, color='green')
    ax.set_xlabel('Time')
    ax.set_ylabel('Activity (a.u.)')
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    return fig
