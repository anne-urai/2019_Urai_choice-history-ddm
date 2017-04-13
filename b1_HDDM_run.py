#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
takes input arguments from stopos
Important: on Cartesius, call module load python/2.7.9 before running
(the only environment where HDDM is installed)
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
# parse input arguments
# ============================================ #

from optparse import OptionParser
usage = "HDDM_run.py [options]"
parser = OptionParser ( usage)
parser.add_option ( "-r", "--run",
        default = 1,
        type = "int",
        help = "Force running the model?" )
parser.add_option ( "-d", "--dataset",
        default = range(0,1),
        type = "int",
        help = "Which dataset, see below" )
parser.add_option ( "-v", "--version",
        default = range(0,8),
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 1,
        type = "int",
        help = "Which trace to run, usually 0-60" )

opts,args       = parser.parse_args()
model_version   = opts.version
d               = opts.dataset
trace_id        = opts.trace_id
runMe           = opts.run

# to avoid errors when plotting on cartesius
# http://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
import matplotlib
matplotlib.use('Agg') # to still plot even when no display is defined
import matplotlib.pyplot as plt

# ============================================ #
# define the function that will do the work
# ============================================ #

def make_model(mypath, model_name, trace_id):

    import os, fnmatch
    import hddm
    import numpy as np

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename

    # get the csv file for this dataset
    filename    = fnmatch.filter(os.listdir(mypath), '*.csv')
    mydata      = hddm.load_csv(os.path.join(mypath, filename[0]))
    mydata      = mydata.dropna(subset=['prevresp']) # dont use trials with nan in prevresp

    # prepare link function for the regression models
    def z_link_func(x, data=mydata):
        return 1 / (1 + np.exp(-(x.values.ravel())))

    def balance_designmatrix(mydata):
        # remove subjects who did not do all conditions
        for i, sj in enumerate(mydata.subj_idx.unique()):
        	sessions = mydata[mydata.subj_idx == sj].session.unique()
        	if len(sessions) < len(mydata.session.unique()):
        		mydata = mydata[mydata.subj_idx != sj] # drop this subject
        return mydata

    # ============================================ #
    # STEP 1. DOES PREVRESP AFFECT DC OR Z?
    # ============================================ #

    if model_name == 'regress_dc_prevresp':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)', 'link_func': z_link_func}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}

        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)', 'link_func': lambda x:x}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)', 'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # STEP 2. SESSION-DEPENDENCE
    # ============================================ #

    if model_name == 'regress_dc_prevresp_sessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp', 'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp*prevrt*C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp*prevrt',
                'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevrt_sessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        if 'transitionprob' in mydata.columns:
            raise ValueError('Do not fit session-specific serial bias on Anke''s data')

        # boundary separation and drift rate will change over sessions
        v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp*prevrt*C(session)',
            'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 +  + prevresp*prevrt*C(transitionprob) + prevresp*prevpupil*C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp*prevrt + prevresp*prevpupil',
                'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevrt_prevpupil_sessions':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        mydata = balance_designmatrix(mydata)

        if 'transitionprob' in mydata.columns:
            raise ValueError('Do not fit session-specific serial bias on Anke''s data')

        # boundary separation and drift rate will change over sessions
        v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp*prevrt*C(session)+ prevpupil*prevrt*C(session)',
            'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    return m

def run_model(m, mypath, model_name, trace_id, nr_samples=15000):

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    m.find_starting_values() # this should help the sampling
    m.sample(nr_samples, burn=nr_samples/3, thin=2, db='pickle',
        dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
    # specify a certain backend? pickle?
    m.save(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)) # save the model to disk

    # ============================================ #
    # save the output values
    # ============================================ #

    results = m.gen_stats() # this seems different from print_stats??
    results.to_csv(os.path.join(mypath, model_name, 'results-md%d.csv'%trace_id))

    # save the DIC for this model
    text_file = open(os.path.join(mypath, model_name, 'DIC-md%d.txt'%trace_id), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

    # ============================================ #
    # save traces
    # ============================================ #

    all_traces = m.get_traces()
    all_traces.to_csv(os.path.join(mypath, model_name, 'all_traces-md%d.csv'%trace_id))

    # ============================================ #
    # plot convergence check
    # ============================================ #

    # plot the traces and posteriors for each parameter
    figpath = os.path.join(mypath, model_name, 'figures-md%d'%trace_id)
    m.plot_posteriors(save=True, path=figpath, format='pdf')

def concat_models(mypath, model_name):

    import os, hddm, time, kabuki
    from IPython import embed as shell

    # ============================================ #
    # APPEND MODELS
    # ============================================ #

    allmodels = []
    print "appending models"
    for trace_id in range(30): # how chains were run?
        model_filename        = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
        modelExists           = os.path.isfile(model_filename)
        print model_filename
        assert modelExists == True # if not, this model has to be rerun
        thism                 = hddm.load(model_filename)
        # now append
        allmodels.append(thism)

    # ============================================ #
    # CHECK CONVERGENCE
    # ============================================ #

    gr = hddm.analyze.gelman_rubin(allmodels)

    # save
    text_file = open(os.path.join(mypath, model_name, 'gelman_rubin.txt'), 'w')
    for p in gr.items():
        text_file.write("%s:%s\n" % p)
        # print a warning when non-convergence is detected
        # Values should be close to 1 and not larger than 1.02 which would indicate convergence problems.
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3731670/
        if abs(p[1]-1) > 0.02:
            print "non-convergence found, %s:%s\n" %p
    text_file.close()
    print "written gelman rubin stats to file"

    # ============================================ #
    # SAVE POINT ESTIMATES
    # ============================================ #

    # now actually concatenate them, see email Gilles
    # THIS ONLY WORKS IF Z HAS BEEN TRANSFORMED!
    m = kabuki.utils.concat_models(allmodels)
    # m.save(os.path.join(mypath, model_name, 'modelfit-combined.model')) # save the model to disk

    results = m.gen_stats() # point estimate for each parameter and subject
    results.to_csv(os.path.join(mypath, model_name, 'results-combined.csv'))

    # save the DIC for this model
    text_file = open(os.path.join(mypath, model_name, 'DIC-combined.txt'), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

    # ============================================ #
    # SAVE TRACES
    # ============================================ #

    # get the names for all nodes that are available here
    group_traces = m.get_group_traces()
    group_traces.to_csv(os.path.join(mypath, model_name, 'group_traces.csv'))

    all_traces = m.get_traces()
    all_traces.to_csv(os.path.join(mypath, model_name, 'all_traces.csv'))

    # plot the traces and posteriors for each parameter
    figpath = os.path.join(mypath, model_name, 'figures-concat')
    m.plot_posteriors(save=True, path=figpath, format='pdf')

# ============================================ #
# PREPARE THE ACTUAL MODEL FITS
# ============================================ #

# which model are we running at the moment?
models = {0: 'regress_dc_prevresp',
    1: 'regress_z_prevresp',
    2: 'regress_dc_z_prevresp',
    3: 'regress_dc_prevresp_sessions',
    4: 'regress_dc_prevresp_prevrt',
    5: 'regress_dc_prevresp_prevrt_sessions',
    6: 'regress_dc_prevresp_prevrt_prevpupil',
    7: 'regress_dc_prevresp_prevrt_prevpupil_sessions'}

datasets = {0: 'RT_RDK', 1: 'MEG', 2: 'Anke_serial'}

# recode
if isinstance(d, int):
    d = range(d,d+1) # makes a list out of an integer
if isinstance(model_version, int):
    model_version = range(model_version, model_version+1)

for dx in d:

    # find path depending on location and dataset
    import os, time
    mypath = os.path.realpath(os.path.expanduser('~/Data/%s/HDDM'%datasets[dx]))

    for vx in model_version:

        # make a folder for the outputs, combine name and time
        thispath = os.path.join(mypath, models[vx])
        if not os.path.exists(thispath):
            os.mkdir(thispath)

        if runMe == True:
            starttime = time.time()
            # get the model specification
            m = make_model(mypath, models[vx], trace_id)
            # now sample and save
            run_model(m, mypath, models[vx], trace_id)
            elapsed = time.time() - starttime
            print( "Elapsed time: %f seconds\n" %elapsed )

        else: # concatenate the different chains
            concat_models(mypath, models[vx])
