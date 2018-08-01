import numpy as np
import matplotlib.pyplot as plt
from IPython import embed as shell

def DDM_traces_get(v=1, z=0.5, dc=0, dc_slope=0, scaling=0.1, stim=0, nr_trials=1000, t=5.0, dt=0.01):
    
    """
    DDM

    z:  Starting point
    v:  Drift rate
    dc: Drift criterion
    ou: Ornstein Uhlenbeck parameter
    """
    
    if stim == 0:
        v = -v
    
    x1 = np.zeros((nr_trials, int(t/dt)))
    x1[:,:] = np.NaN
    x1[:,0] = z
    for i in range((int(t/dt))-1):
        x1[:,i+1] = x1[:,i] + ((v + dc + (dc_slope * dt * i) ) * dt) + (np.random.normal(0,1,nr_trials)*np.sqrt(dt)*scaling)
    return x1
    
def OU_traces_get(v, ou, dc, z, pre_generated=False, stim=0, nr_trials=1000, timesteps=200,):
    
    """
    LCA
    nonlinear == True --> LCA
    nonlinear == False --> OU
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
            x1[:,i+1] = x1[:,i] + v[0] + dc[0] + (ou[0]*x1[:,i]) + np.random.normal(0,0.1,nr_trials)
            x2[:,i+1] = x2[:,i] + v[1] + dc[1] + (ou[1]*x2[:,i]) + np.random.normal(0,0.1,nr_trials)
    return x1, x2

def LCA_traces_get(v, k, w, dc, z, pre_generated=False, stim=0, nr_trials=1000, timesteps=200,):
    
    """
    LCA
    nonlinear == True --> LCA
    nonlinear == False --> OU
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
            x1[:,i+1] = np.clip(x1[:,i] + v[0] + dc[0] - (k[0]*x1[:,i]) - (w[1]*x2[:,i]) + np.random.normal(0,0.1,nr_trials), a_min=0, a_max=1e6)
            x2[:,i+1] = np.clip(x2[:,i] + v[1] + dc[1] - (k[1]*x2[:,i]) - (w[0]*x1[:,i]) + np.random.normal(0,0.1,nr_trials), a_min=0, a_max=1e6)
                    
    return x1, x2

def one_accumulater_traces_apply_bounds(x1, a=0.15, b0_collapse=0, b1_collapse=0, dt=0.01):
    
    rt = np.zeros(x1.shape[0])
    response = np.zeros(x1.shape[0])
    for i in range(x1.shape[0]):
        bound1 = a
        bound0 = 0
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
