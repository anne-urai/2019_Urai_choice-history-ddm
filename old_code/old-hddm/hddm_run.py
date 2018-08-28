#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
adapted from JW de Gee
"""

from joblib import Parallel, delayed
nr_traces = 3
model_name = 'basic_stimcoding'

def run_model(trace_id):

    import os, pickle, sys, time, hddm
    model_name  = 'basic_stimcoding'
    nsamples    = 10000 # 50.000 for real results?

    # find path depending on local/klimag
    usr = os.environ.get('USER')
    if usr in ['anne']:
        mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
    if usr in ['aurai']:
        mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'

     # make a folder for the outputs, combine name and time
     # currentTime = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
     # thispath    = os.path.join(mypath, model_name + '_' + currentTime)
    thispath    = os.path.join(mypath, model_name)
    if not os.path.exists(thispath):
        os.mkdir(thispath)

    # ============================================ #
    # load in data
    # ============================================ #

    mydata = hddm.load_csv(os.path.join(mypath, '2ifc_data_hddm.csv'))

    # specify the model
    m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
        drift_criterion=True, bias=True,
        include=('sv'), group_only_nodes=['sv'],
        depends_on={'t':['sessionnr'], 'v':['sessionnr'],
        'a':['sessionnr'], 'dc':['sessionnr'], 'z':['sessionnr']},
        p_outlier=.05)

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    model_filename = os.path.join(thispath, 'modelfit-md%d'%trace_id)
    m.sample(nsamples, burn=nsamples/4, thin=3, db='pickle',
        dbname=model_filename)
    m.save(model_filename) # save the model, see if this works with pickle

    # ============================================ #
    # save the output values
    # ============================================ #

    results = m.gen_stats()
    results.to_csv(os.path.join(thispath, 'results-md%d.csv'%trace_id))

    # save the DIC for this model
    text_file = open(os.path.join(thispath, 'DIC-md%d.txt'%trace_id), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

    # plot some output stuff in figures subfolder
    figpath = os.path.join(thispath, 'figures-md%d'%trace_id)
    if not os.path.exists(figpath):
        os.mkdir(figpath)

    m.plot_posteriors(save=True, path=figpath, format='pdf')
    # m.plot_posterior_predictive(samples=10, bins=100, num_subjs=10,
    # save=False, path=figpath, format='pdf')

# ============================================ #
# and go!
# ============================================ #

Parallel(n_jobs=nr_traces, verbose=10)(delayed(run_model)(k) for k in range(nr_traces))
