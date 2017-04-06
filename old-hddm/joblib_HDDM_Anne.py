#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
adapted from JW de Gee
uses joblib instead of parallel python
"""

# ============================================ #
# HDDM cheat sheet
# ============================================ #

# v     = drift rate
# a     = boundary separation
# t     = nondecision time
# z     = starting point
# dc    = drift driterion
# sv    = inter-trial variability in drift-rate
# st    = inter-trial variability in non-decision time
# sz    = inter-trial variability in starting-point

# ============================================ #
# define the function that will do the work
# ============================================ #

def run_model(model_name, trace_id, nr_samples=10000):

    import hddm, os
    import matplotlib.pyplot as plt

    # find path depending on local/klimag
    usr = os.environ.get('USER')
    if usr in ['anne']:
        mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
    if usr in ['aurai']:
        mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'

    # make a folder for the outputs, combine name and time
    thispath = os.path.join(mypath, model_name)
    if not os.path.exists(thispath):
        os.mkdir(thispath)

    # get the csv
    mydata = hddm.load_csv(os.path.join(mypath, '2ifc_data_hddm.csv'))

    # specify the model
    m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
        drift_criterion=True, bias=True,
        include=('sv'), group_only_nodes=['sv'],
        depends_on={'t':['sessionnr'], 'v':['sessionnr'],
        'a':['sessionnr'], 'dc':['sessionnr', 'prevresp'], 'z':['sessionnr']},
        p_outlier=.05)

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    model_filename = os.path.join(thispath, 'modelfit-md%d'%trace_id)
    m.sample(nr_samples, burn=nr_samples/4, thin=3, db='pickle',
        dbname=model_filename)
    m.save(model_filename) # save the model as pickle

    # ============================================ #
    # save the output values
    # ============================================ #

    results = m.gen_stats() # this seems different from print_stats??
    results.to_csv(os.path.join(thispath, 'results-md%d.csv'%trace_id))

    # save the DIC for this model
    text_file = open(os.path.join(thispath, 'DIC-md%d.txt'%trace_id), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

    return m

# ============================================ #
# set up this model
# ============================================ #

model_name      = 'prevresp_dc_stimcoding'
nr_samples      = 50 # 50.000 for real results?
useParallel     = True # use parallel processing?
nr_traces       = 3

# ============================================ #
# load in data
# ============================================ #

import os, sys, time, hddm
if useParallel:
    from joblib import Parallel, delayed

# find path depending on local/klimag
usr = os.environ.get('USER')
if usr in ['anne']:
    mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
if usr in ['aurai']:
    mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'
thispath      = os.path.join(mypath, model_name)

# ============================================ #
# and go!
# ============================================ #

start_time  = time.time()

if useParallel: # run in parallel
    models = Parallel(n_jobs=nr_traces, verbose=100) \
        (delayed(run_model)(model_name, trace_id, nr_samples) \
        for trace_id in range(nr_traces))
else:
    models = []
    for trace_id in range(nr_traces): # run the models serially
        thism = run_model(model_name, trace_id)
        models.append(thism)

# ============================================ #
# post-process
# ============================================ #

# for trace_id in range(nr_traces):
#     # plot some output stuff in figures subfolder
#     figpath = os.path.join(thispath, 'figures-md%d'%trace_id)
#     if not os.path.exists(figpath):
#         os.mkdir(figpath)
#
#     models[trace_id].plot_posteriors(save=True, path=figpath, format='pdf')
#     plt.close('all') # this will leave figures open, make sure to close them all
#     # m.plot_posterior_predictive(samples=10, bins=100, num_subjs=10,
#     # save=False, path=figpath, format='pdf')

# gelman rubic, only makes sense when several models were run
gr = hddm.analyze.gelman_rubin(models)
text_file = open(os.path.join(thispath, 'gelman_rubic.txt'), 'w')
for p in gr.items():
     text_file.write("%s:%s\n" % p)
text_file.close()

print "DONE! Time elapsed: ", time.time() - start_time, "s"
