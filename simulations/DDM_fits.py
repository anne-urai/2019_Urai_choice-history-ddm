#!/usr/bin/env python
# encoding: utf-8
"""
Created by Jan Willem de Gee on 2011-02-16.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""
import os, sys, pickle, time
import datetime
import collections
import math
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import hddm
import kabuki
from IPython import embed as shell
from joblib import Parallel, delayed

matplotlib.rcParams['pdf.fonttype'] = 42
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

data_dir = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/fits/')
model_dir = os.path.expanduser('~/projects/2018_Urai_choice-history-ddm/fits/')

global parallel
parallel = False

datasets = [
            "2018_ddm_data_1",                       # 0
            "2018_ddm_data_2",                       # 1
            "2018_ddm_data_3",                       # 2
            "2018_ddm_data_4",                       # 3
            "2018_ddm_data_5",                       # 4
            "2018_ou_data_1",                        # 5
            "2018_ou_data_2",                        # 6
            "2018_ou_data_3",                        # 7
            "2018_ddm_autocorr_data_1",                       # 8
            "2018_ddm_autocorr_data_2",                       # 9
            "2018_ddm_autocorr_data_3",                       # 10
            ]

datasets = ["2018_ddm_autocorr_data_1",                       # 8
            "2018_ddm_autocorr_data_2",                       # 9
            "2018_ddm_autocorr_data_3"]

def aic(self):
    k = len(self.get_stochastics())
    logp = sum([x.logp for x in self.get_observeds()['node']])  
    return 2 * k - 2 * logp

def bic(self):
    k = len(self.get_stochastics())
    n = len(self.data)
    logp = sum([x.logp for x in self.get_observeds()['node']])
    return -2 * logp + k * np.log(n)
    
def fit_ddm_per_group(data, model, model_dir, model_name, samples=5000, burn=1000, thin=1, n_models=3, n_jobs=12):
    
    if parallel:
        res = Parallel(n_jobs=n_jobs)(delayed(fit_ddm_hierarchical)(df, model, model_dir, model_name, samples, burn, thin, model_id) for model_id in range(n_models))
    else:
        fit_ddm_hierarchical(df, model, model_dir, model_name, samples, burn, thin, 0)

def fit_ddm_hierarchical(data, model, model_dir, model_name, samples=5000, burn=1000, thin=1, model_id=0):
    
    exec('global m; m = {}'.format(model))
    m.find_starting_values()
    m.sample(samples, burn=burn, thin=thin, dbname=os.path.join(model_dir, '{}_{}.db'.format(model_name, model_id)), db='pickle')
    m.save(os.path.join(model_dir, '{}_{}.hddm'.format(model_name, model_id)))

    # GET INDIVIDUAL SUBJECT PARAMETERS, FLAT FORMAT
    print("saving stats")
    results = m.gen_stats() # point estimate for each parameter and subject
    results.to_csv(os.path.join(model_dir, '{}_{}_results.csv'.format(model_name, model_id)))

    # FIT MODEL COMPARISON INDICES
    df = dict()
    df['dic_original'] = [m.dic]
    df['aic'] = [aic(m)]
    df['bic'] = [bic(m)]
    df2 = pd.DataFrame(df)
    df2.to_csv(os.path.join(model_dir, '{}_{}_model_comparison.csv'.format(model_name, model_id)))

    return m

def load_ddm_per_group(model_dir, model_name, n_models=3):
    
    models = [kabuki.utils.load(os.path.join(model_dir, '{}_{}.hddm'.format(model_name, model_id))) for model_id in range(n_models)] 
    
    return models

def fit_ddm_per_subject(data, model, model_dir, model_name, n_runs=5, n_jobs=12):
    
    # res = Parallel(n_jobs=n_jobs, backend='loky')(delayed(fit_ddm_subject)(subj_data, subj_idx, model, model_dir, model_name, n_runs) for subj_idx, subj_data in df.groupby('subj_idx'))
    res = []
    for subj_idx, subj_data in df.groupby('subj_idx'):
        res.append(fit_ddm_subject(subj_data, subj_idx, model, model_dir, model_name, n_runs))
        
    res = pd.concat(res, axis=0)
    res.to_csv(os.path.join(model_dir, '{}_params_flat.csv'.format(model_name)))
    return res

def fit_ddm_subject(data, subj_idx, model, model_dir, model_name, n_runs=5):

    data = data.loc[data["subj_idx"]==subj_idx,:]
    exec('global m; m = {}'.format(model))

    # m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
    # p_outlier=0.05)
    
    # optimize:
    # m.approximate_map()
    m.optimize('gsquare', quantiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], n_runs=n_runs)
    res = pd.concat((pd.DataFrame([m.values], index=[subj_idx]), pd.DataFrame([m.bic_info], index=[subj_idx])), axis=1)

    try:
        res['aic'] = m.aic
        res['bic'] = m.bic
    except:
        res['aic'] = np.nan
        res['bic'] = np.nan

    # compare with bic_info from kabuki
    bic_info = m.bic_info
    res['bic_info'] = bic_info['bic']

    return res

def load_ddm_per_subject(model_dir, model_name):
    
    return pd.read_csv(os.path.join(model_dir, '{}_params_flat.csv'.format(model_name))).drop('Unnamed: 0', 1)

# dataset = 4
# version = 0

run = True
for ds, dataset in enumerate(datasets):
    # for version in [0,1,2,3]:
    for version in [4,5,6]:
        
        # load data:
        print(os.path.join(data_dir, '{}.csv'.format(datasets[ds])))
        df = pd.read_csv(os.path.join(data_dir, '{}.csv'.format(datasets[ds])))

        # stimcoding?
        stimcoding = True
        df.rename(columns={'choice_a':'response'}, inplace=True)

        # variables:
        subjects = np.unique(df.subj_idx)
        nr_subjects = len(subjects)

        # model:
        if version == 0:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=False, bias=False, p_outlier=0)"
        if version == 1:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=False, bias=True, p_outlier=0)"
        elif version == 2:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=False, p_outlier=0)"
        elif version == 3:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0)"

        elif version == 4:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0, depends_on={'z':['prevresp']})"
        elif version == 5:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0, depends_on={'dc':['prevresp']})"
        elif version == 6:
            model = "hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0, " \
                    "depends_on={'dc':['prevresp'], 'z':['prevresp']})"

        # model_name:
        model_name = '{}_{}'.format(datasets[ds], version)
    
        # create figure dir:
        fig_dir = os.path.join(model_dir, model_name)
        
        # fit model:
        if run:
            
            print("fitting {}".format(model_name))
            n_jobs = 24

            # flat:
            results = fit_ddm_per_subject(df, model, model_dir, model_name, n_runs=5, n_jobs=n_jobs)

            # # hierarchical:
            # results = fit_ddm_per_group(df, model, model_dir, model_name, n_jobs=2)

            # save:
            # results.to_csv(os.path.join(fig_dir, 'results.csv'))

            # # barplot:
            # fig = plt.figure(figsize=(6,2))
            # ax = fig.add_subplot(111)
            # sns.barplot(data=results.loc[:,['a', 'v', 't', 'z', 'dc']], ax=ax)
            # sns.despine(offset=5, trim=True)
            # plt.tight_layout()
            # fig.savefig(os.path.join(fig_dir, 'bars.pdf'))
            
            # # barplot only z and dc:
            # # results['z'] = results['z'] - 0.5
            # fig = plt.figure(figsize=(1.25,1.6))
            # ax = fig.add_subplot(111)
            # sns.barplot(data=results.loc[:,['z', 'dc']], ax=ax)
            # ax.set_ylim(0,0.8)
            # sns.despine(offset=5, trim=True)
            # plt.tight_layout()
            # fig.savefig(os.path.join(fig_dir, 'bars2.pdf'))