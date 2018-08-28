#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
adapted from JW de Gee
uses IPython parallel instead of parallel python
"""

import hddm, sys, time
from datetime import datetime
import ipyparallel as ipp # dont use parallel python but run within ipython instead
from IPython.display import clear_output

useIPP      = True # use parallel processing?
nr_traces   = 3

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

def run_model(trace_id):

    import os, pickle, sys, time, hddm
    model_name  = 'basic_stimcoding'
    nsamples    = 10 # 50.000 for real results?

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
    # figpath = os.path.join(thispath, 'figures-md%d'%trace_id)
    # if not os.path.exists(figpath):
    #     os.mkdir(figpath)

    # m.plot_posteriors(save=True, path=figpath, format='pdf')
    # m.plot_posterior_predictive(samples=10, bins=100, num_subjs=10,
    # save=False, path=figpath, format='pdf')

    # print "model fit completed"
    return m

# ============================================ #
# function to track parallel output
# from https://groups.google.com/forum/#!topic/hddm-users/lsmRwwB_9wY
# ============================================ #

def wait_watching_stdout(ar, dt=5):
    ## ar: vmap output of the models being run
    ## dt: number of seconds between checking output
    while not ar.ready():
        stdouts = ar.stdout
        if not any(stdouts):
            continue
        # clear_output doesn't do much in terminal environments
        clear_output()
        print ""
        for out in ar.stdout: print(out);
        sys.stdout.flush()
        time.sleep(dt)

# ============================================ #
# and go!
# ============================================ #

if useIPP: # run in parallel

    print "starting parallel client"

    c       = ipp.Client() # start the client
    rview   = c[range(nr_traces)] # use one core for each instantiation of the model
    # rview   = c.load_balanced_view()
    # only map_async allows tracking of the output

    print "mapping jobs"
    jobs    = rview.map_async(run_model, range(nr_traces))
    # wait_watching_stdout(jobs) # print output to command line to track progress
    print "waiting for jobs to complete"

    models  = jobs.get() # will wait for the output
    print "parallel running succeeded"

else:
    models = []
    for t in range(nr_traces): # run the models serially
        thism = run_model(t)
        models.append(thism)

# from IPython import embed as shell
# shell()
# gelman rubic, only makes sense when several models were run
gr = hddm.analyze.gelman_rubin(models)
text_file = open(os.path.join(thispath, 'gelman_rubic.txt'), 'w')
for p in gr.items():
     text_file.write("%s:%s\n" % p)
text_file.close()

print "DONE! Time elapsed: ", time.time() - start_time, "s"
