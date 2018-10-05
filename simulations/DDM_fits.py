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

import json

from jw_tools import myfuncs
from jw_tools import ddm_tools

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

base_dir = os.path.expanduser('~/Desktop/simulations/')
data_dir = os.path.join(base_dir, 'ddm_fits_data')
model_dir = os.path.join(base_dir, 'ddm_fits_model')

datasets = [
            "2018_ddm_data_1",                      # 0
            "2018_ddm_data_2",                      # 1
            "2018_ddm_data_3",                      # 2
            "2018_ou_data_1",                       # 3
            "2018_ou_data_2",                       # 4
            "2018_ou_data_3",                       # 5
            "2018_ou_data_4",                       # 6
            "2018_ou_data_5",                       # 7
            "2018_ou_data_6",                       # 8
            "2018_lca_data_1",                      # 9
            "2018_lca_data_2",                      # 10
            "2018_lca_data_3",                      # 11
            "2018_lca_data_4",                      # 12
            ]
    
def fit_ddm_per_group(data, model, model_dir, model_name, samples=5000, burn=1000, thin=1, n_models=3, n_jobs=12):
    
    res = Parallel(n_jobs=n_jobs)(delayed(fit_ddm_hierarchical)(df, model, model_dir, model_name, samples, burn, thin, model_id) for model_id in range(n_models))

def fit_ddm_hierarchical(data, model, model_dir, model_name, samples=5000, burn=1000, thin=1, model_id=0):
    
    exec('global m; m = {}'.format(model))
    m.find_starting_values()
    m.sample(samples, burn=burn, thin=thin, dbname=os.path.join(model_dir, '{}_{}.db'.format(model_name, model_id)), db='pickle')
    m.save(os.path.join(model_dir, '{}_{}.hddm'.format(model_name, model_id)))
    
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
    
    import hddm

    data = data.loc[data["subj_idx"]==subj_idx,:]
    exec('global m; m = {}'.format(model))

    # optimize:
    m.optimize('gsquare', quantiles=analysis_info['quantiles'], n_runs=n_runs)
    res = pd.concat((pd.DataFrame([m.values], index=[subj_idx]), pd.DataFrame([m.bic_info], index=[subj_idx])), axis=1)
    
    return res

def load_ddm_per_subject(model_dir, model_name):
    
    return pd.read_csv(os.path.join(model_dir, '{}_params_flat.csv'.format(model_name))).drop('Unnamed: 0', 1)

# dataset = 4
# version = 0

run = True
for ds in [0,1,2,3,4,5,6,7,8,9,10,11,12]:
    for version in [0]:
        
        # load analysis info:
        with open(os.path.join(data_dir, '{}.json'.format(datasets[ds]))) as json_data:
            analysis_info = json.load(json_data)
        
        # load data:
        try:
            data = pd.read_csv(os.path.join(data_dir, analysis_info['data_file'])).drop('Unnamed: 0', 1)
        except:
            data = pd.read_csv(os.path.join(data_dir, analysis_info['data_file']))
        
        # stimcoding?
        if analysis_info['stimcoding'] == "True":
            stimcoding = True
            data.rename(columns={'choice_a':'response'}, inplace=True)
        else:
            stimcoding = False
            data.rename(columns={'correct':'response'}, inplace=True)
        
        # variables:
        subjects = np.unique(data.subj_idx)
        nr_subjects = len(subjects)
        model = analysis_info["model"][version]
        
        # run:
        for split_by in analysis_info['split_by']:
            
            # model_name:
            model_name = '{}_{}_{}'.format(analysis_info["model_name"], split_by, version)
        
            # create figure dir:
            fig_dir = model_dir = os.path.join(base_dir, 'ddm_fits_figs', model_name)
            try:
                os.system('mkdir {}'.format(fig_dir))
                os.system('mkdir {}'.format(os.path.join(fig_dir, 'diagnostics')))
            except:
                pass
            
            # prepare dataframe:
            df = data.copy()
                        
            # fit model:
            if run:
                
                print("fitting {}".format(model_name))
                n_jobs = 4
                
                # # hierarchical:
                # results = fit_ddm_per_group(data, model, model_dir, model_name, samples=5000, burn=1000, thin=1, n_models=3, n_jobs=n_jobs)
                
                # flat:
                results = fit_ddm_per_subject(df, model, model_dir, model_name, n_runs=5, n_jobs=n_jobs)
                
            # for fit_type in ['flat', 'hierarchical']:
            for fit_type in ['flat']:
            
                if fit_type == 'flat':
                    results = load_ddm_per_subject(model_dir, model_name)
                else:
                    models = load_ddm_per_group(model_dir, model_name, n_models=3)
                    
                    # gelman rubic:
                    gr = hddm.analyze.gelman_rubin(models)
                    text_file = open(os.path.join(fig_dir, 'diagnostics', 'gelman_rubic_{}.txt'.format(fit_type)), 'w')
                    for p in gr.items():
                        text_file.write("%s:%s\n" % p)
                    text_file.close()

                    # dic:
                    text_file = open(os.path.join(fig_dir, 'diagnostics', 'DIC_{}.txt'.format(fit_type)), 'w')
                    for i, m in enumerate(models):
                        text_file.write("Model {}: {}\n".format(i, m.dic))
                    text_file.close()
                    
                    # posteriors:
                    m.plot_posteriors(save=True, path=os.path.join(fig_dir, 'diagnostics'), format='pdf')
                    
                    # dataframe:
                    m = models[1]
                    results = m.gen_stats()['50q'].reset_index()
                    results = results.loc[['subj' in c for c in results["index"]],:]
                    a = np.array([results.iloc[i]["index"].split('.')[0] for i in range(results.shape[0])])
                    _, idx = np.unique(a, return_index=True)
                    cols = a[np.sort(idx)]
                    cols = np.array([c.replace('_subj', '') for c in cols])
                    results = pd.DataFrame(np.vstack([np.array(results.loc[np.array([str(subj_idx) == c.split('.')[-1] for c in results["index"]]),:]["50q"]) for subj_idx in subjects]))
                    results.columns = cols
                try:
                    params = results.drop(['bic', 'likelihood', 'penalty'], 1)
                except:
                    params = results.copy()
                params = params.loc[:,~np.array(['z_trans' in c for c in params.columns])]
                params.to_csv(os.path.join(fig_dir, 'params_{}.csv'.format(fit_type)))
                
                # barplot:
                fig = plt.figure(figsize=(6,2))
                ax = fig.add_subplot(111)
                sns.barplot(data=params, ax=ax)
                sns.despine(offset=5, trim=True)
                plt.tight_layout()
                fig.savefig(os.path.join(fig_dir, 'bars_{}.pdf'.format(fit_type)))
                
                # barplot only z and dc:
                params['z'] = params['z'] - 0.5
                fig = plt.figure(figsize=(1.25,1.6))
                ax = fig.add_subplot(111)
                sns.barplot(data=params.loc[:,['z', 'dc']], ax=ax)
                ax.set_ylim(0,0.8)
                sns.despine(offset=5, trim=True)
                plt.tight_layout()
                fig.savefig(os.path.join(fig_dir, 'bars2_{}.pdf'.format(fit_type)))