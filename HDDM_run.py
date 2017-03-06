#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
takes input arguments from stopos
Important: on Cartesius, call module load python2.7.9 before running
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
parser.add_option ( "-d", "--dataset",
        default = 0,
        type = "int",
        help = "Which dataset, see below" )
parser.add_option ( "-v", "--version",
        default = 0,
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 1,
        type = "int",
        help = "Which trace to run, usually 1-3" )

opts,args       = parser.parse_args()
model_version   = opts.version
d               = opts.dataset
trace_id        = opts.trace_id

# ============================================ #
# define the function that will do the work
# ============================================ #

def run_model(mypath, model_name, trace_id, nr_samples=50):

    import os, fnmatch
    import hddm

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename

    # get the csv file for this dataset
    filename = fnmatch.filter(os.listdir(mypath), '*.csv')
    mydata = hddm.load_csv(os.path.join(mypath, filename[0]))
    print mydata.head(n=10) # show the data

    # specify the model
    if model_name == 'stimcoding':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05) #,
        #    depends_on={'v':['session']})

    elif model_name == 'stimcoding_prevresp_dc':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            depends_on={'v':['session'], 'dc':['prevresp']})

    elif model_name == 'stimcoding_prevresp_z':
        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            depends_on={'v':['session'], 'z':['prevresp']})

    elif model_name == 'regress_dc':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        v_reg = {'model': 'v ~ 1 + stimulus + stimulus:session', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg, include='z', group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp']) # dont use trials with nan in prevresp

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        v_reg = {'model': 'v ~ 1 + stimulus + stimulus:session + prevresp', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg, include='z', group_only_regressors=False, p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevpupil_prevrt':
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1         # recode the stimuli into signed
        mydata = mydata.dropna(subset=['prevresp', 'prevpupil']) # dont use trials with nan in prevresp or prevpupil

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        v_reg = {'model': 'v ~ 1 + stimulus + stimulus:session + prevresp + prevresp:prevrt + prevresp:prevpupil', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg, include='z', group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    m.sample(nr_samples, burn=nr_samples/10, thin=2, db='pickle',
        dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
    m.print_stats() # just for display in command window
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

# ============================================ #
# run one model per job
# ============================================ #

# which model are we running at the moment?
models = {0: 'stimcoding',
    1: 'stimcoding_prevresp_dc',
    2: 'stimcoding_prevresp_z',
    3: 'regress_dc',
    4: 'regress_dc_prevresp',
    5: 'regress_dc_prevresp_prevpupil_prevrt'}

datasets = {0: 'RT_RDK',
    1: 'MEG-PL'}

# find path depending on location and dataset
import os, time
mypath = os.path.expanduser('~/Data/%s/HDDM'%datasets[d])

# make a folder for the outputs, combine name and time
thispath = os.path.join(mypath, models[model_version])
time.sleep(trace_id) # make sure this doesn't crash the script if multiple instances of the model are submitted simultaneously
if not os.path.exists(thispath):
    os.mkdir(thispath)

# and... go!
starttime = time.time()
run_model(mypath, models[model_version], trace_id)
elapsed = time.time() - starttime
print( "Elapsed time: %f seconds\n" %elapsed )

# and plot
# import HDDM_plotOutput
# HDDM_plotOutput.plot_model(mypath, models[model_version], trace_id)
