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
# define the function that will do the workq
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
    # print mydata.head(n=5) # show the data

    # prepare link functions for the regression models
    def z_link_func(x, data=mydata):
        return 1 / (1 + np.exp(-(x.values.ravel())))

    # ============================================ #
    # specify the model
    # ============================================ #

    if model_name == 'stimcoding':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'v':['session']})

    elif model_name == 'stimcoding_prevresp_dc':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'v':['session'], 'dc':['prevresp']})

    elif model_name == 'stimcoding_prevresp_z':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'v':['session'], 'z':['prevresp']})

    elif model_name == 'stimcoding_prevresp_dc_z':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'v':['session'], 'z':['prevresp'], 'dc':['prevresp']})

    elif model_name == 'regress_dc':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed

        v_reg = {'model': 'v ~ 1 + stimulus*session', 'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp']) # dont use trials with nan in prevresp

        v_reg = {'model': 'v ~ 1 + stimulus*session + prevresp', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevpupil_prevrt':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp', 'prevpupil']) # dont use trials with nan in prevresp or prevpupil

        v_reg = {'model': 'v ~ 1 + stimulus*session + prevresp + prevpupil:prevresp + prevrt:prevresp', 'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_z_prevresp':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp']) # dont use trials with nan in prevresp

        z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus*session', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_z_prevresp_prevpupil_prevrt':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp', 'prevpupil']) # dont use trials with nan in prevresp or prevpupil

        z_reg = {'model': 'z ~ 1 + prevresp + prevpupil:prevresp + prevrt:prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus*session', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp']) # dont use trials with nan in prevresp or prevpupil

        z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus*session + prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp_prevpupil_prevrt':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp', 'prevpupil']) # dont use trials with nan in prevresp or prevpupil

        z_reg = {'model': 'z ~ 1 + prevresp + prevpupil:prevresp + prevrt:prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus*session + prevresp + prevpupil:prevresp + prevrt:prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    return m

def run_model(m, mypath, model_name, trace_id, nr_samples=10000):

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    m.sample(nr_samples, burn=nr_samples/3, thin=2)

    #, db='pickle',
    #    dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
    # m.print_stats() # just for display in command window
    # specify a certain backend? pickle?
    # m.save(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)) # save the model to disk

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
    # MANUALLY APPEND MODELS
    # ============================================ #

    # shell()

    # if this was a regression model with a custom link function for z
    if model_name.find('regress_z') > -1:

        # # push z through the reverse logistic function to make sense
        # def z_reverse_link_func(x):
        #     import numpy as np
        #     return 1 / (1 + np.exp(-(x)))
        #     # return np.exp(x) / (1 + np.exp((x)))

        allmodels = []
        for trace_id in range(60): # 15 models were run

            # manually this specific database (which has been pickled)
            db_filename = os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id)
            thism       = make_model(mypath, model_name, trace_id) # first, get the stuff we want
            thism.load_db(db_filename, db='pickle')

            # a = thism.get_stochastics() # pandas dataframe
            # rownames = list(a.index)
            # for idx, znode in enumerate(thism.get_stochastics().node):
            #     # determine whether this node has z
            #     if 'z' in rownames[idx]:
            #         print(znode)
            #         # print(rownames[idx])
            #         znode.trace._trace[0] = z_reverse_link_func(znode.trace._trace[0])
            #         thism.get_stochastics().node[idx] = znode

            # then append as usual
            allmodels.append(thism)

    # ============================================ #
    # APPEND MODELS
    # ============================================ #

    else:
        allmodels = []
        print "appending models"
        for trace_id in range(60): # 15 models were run
            model_filename              = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
            modelExists                 = os.path.isfile(model_filename)
            print model_filename
            assert modelExists == True # if not, this model has to be rerun
            thism                       = hddm.load(model_filename)
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

    if model_name.find('regress_z') < 0:

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

# ============================================ #
# PREPARE THE ACTUAL MODEL FITS
# ============================================ #

# which model are we running at the moment?
models = {0: 'stimcoding_prevresp_dc',
    1: 'stimcoding_prevresp_z',
    2: 'stimcoding_prevresp_dc_z',
    3: 'regress_dc_prevresp',
    4: 'regress_dc_prevresp_prevpupil_prevrt',
    5: 'regress_z_prevresp',
    6: 'regress_z_prevresp_prevpupil_prevrt',
    7: 'regress_dc_z_prevresp',
    8: 'regress_dc_z_prevresp_prevpupil_prevrt'}

datasets = {0: 'RT_RDK', 1: 'MEG-PL', 2: 'MEG-PL-S1', 3: 'MEG-PL-S2'}

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
            m = make_model(mypath, models[vx], trace_id) # first, get the stuff we want
            # now sample and save
            run_model(m, mypath, models[vx], trace_id)
            elapsed = time.time() - starttime
            print( "Elapsed time: %f seconds\n" %elapsed )

        else: # concatenate the different chains
            concat_models(mypath, models[vx])
