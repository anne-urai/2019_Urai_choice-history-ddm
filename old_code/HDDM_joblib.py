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

def run_model(mypath, model_version, trace_id, nr_samples=50000):

    if model_version == 1:
        model_name = 'stimcoding'
    elif model_version == 2:
        model_name = 'prevresp_z'
    elif model_version == 3:
        model_name = 'prevresp_dc'
    elif model_version == 4:
        model_name = 'prevresp_prevrt_z'
    elif model_version == 5:
        model_name = 'prevresp_prevrt_dc'
    elif model_version == 6:
        model_name = 'prevresp_prevpupil_z'
    elif model_version == 7:
        model_name = 'prevresp_prevpupil_dc'

    import os
    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    modelExists     = os.path.isfile(model_filename)

    if modelExists:
        print "model already exists, skipping"
    else:
        import hddm
        # get the csv
        mydata = hddm.load_csv(os.path.join(mypath, '2ifc_data_hddm.csv'))

        # specify the model
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'t':['sessionnr'], 'v':['sessionnr'],
            'a':['sessionnr'], 'dc':['sessionnr'], 'z':['sessionnr', 'prevresp']},
            p_outlier=.05)

        # ============================================ #
        # do the actual sampling
        # ============================================ #

        m.sample(nr_samples, burn=nr_samples/4, thin=3, db='pickle',
            dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
        m.save(model_filename) # save the model to disk

        # ============================================ #
        # save the output values
        # ============================================ #

        results = m.gen_stats() # this seems different from print_stats??
        results.to_csv(os.path.join(mypath, model_name, 'results-md%d.csv'%trace_id))

        # save the DIC for this model
        text_file = open(os.path.join(mypath, model_name, 'DIC-md%d.txt'%trace_id), 'w')
        text_file.write("Model {}: {}\n".format(trace_id, m.dic))
        text_file.close()

    # dont return model object, can't be pickled so Parallel will hang
    return trace_id

# ============================================ #
# set up this model
# ============================================ #

model_name      = 'prevresp_z_stimcoding'
nr_samples      = 50000 # 50.000 for real results?
nr_traces       = 3

# find path depending on local/klimag
import os

usr = os.environ.get('USER')
if usr in ['anne']:
    mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
if usr in ['aurai']:
    mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'

# make a folder for the outputs, combine name and time
thispath = os.path.join(mypath, model_name)
if not os.path.exists(thispath):
    os.mkdir(thispath)

# ============================================ #
# run models in parallel
# ============================================ #

from joblib import Parallel, delayed

Parallel(n_jobs=nr_traces, verbose=10) \
    (delayed(run_model)(mypath, model_name, trace_id, nr_samples) \
    for trace_id in range(nr_traces))

# ============================================ #
# post-processing
# ============================================ #

import hddm
import matplotlib.pyplot as plt

print "HDDM imported, starting post-processing"
models = []
for trace_id in range(nr_traces): # run the models serially
    thism = hddm.load(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id))
    print os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)

    # plot some output stuff in figures subfolder
    figpath = os.path.join(mypath, model_name, 'figures-md%d'%trace_id)
    if not os.path.exists(figpath):
        os.mkdir(figpath)
    thism.plot_posteriors(save=True, path=figpath, format='pdf')
    plt.close('all') # this will leave figures open, make sure to close them all
    models.append(thism)

# gelman rubic on the list of models
gr = hddm.analyze.gelman_rubin(models)
text_file = open(os.path.join(mypath, model_name, 'gelman_rubic.txt'), 'w')
for p in gr.items():
     text_file.write("%s:%s\n" % p)
text_file.close()
