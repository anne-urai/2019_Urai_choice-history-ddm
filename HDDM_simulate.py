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
usage = "runHDDM_stopos.py [options]"
parser = OptionParser ( usage)
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
trace_id        = opts.trace_id

# ============================================ #
# define the function that will do the work
# ============================================ #

def run_model(mypath, model_name, trace_id, nr_samples=1000):

    import os
    import hddm

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename
    modelExists     = os.path.isfile(model_filename)

    #if not modelExists:

    # get the csv
    mydata = hddm.load_csv(os.path.join(mypath, 'rtrdk_data_allsj.csv'))

    # specify the model
    if model_name == 'stimcoding':

        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05)
            # depends_on={'v':['sessionnr'], 'dc':['sessionnr']},

    elif model_name == 'regress_dc':

        # recode the stimuli into signed
        mydata.ix[mydata['stimulus']==0,'stimulus'] = -1

        # question: should response have the choice or accuracy coding?

        # this only works if all subjects have data in all conditions
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        # include a term for serial bias and its modulation by pupil/rt
        # v_reg = {'model': 'v ~ 1 + C(sessionnr > 2):stimulus + prevresp + prevresp*prevpupil + prevresp*prevrt',
        #    'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg, include='z', group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    m.sample(nr_samples, burn=nr_samples/10, thin=1, db='pickle',
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

# ============================================ #
# run one model per job
# ============================================ #

# find path depending on location
import os, time
mypath = '/Users/anne/Data/RT_RDK/HDDM'

# which model are we running at the moment?
models = {0: 'stimcoding',
    1: 'regress_dc',
    2: 'regress_z'}

# make a folder for the outputs, combine name and time
thispath = os.path.join(mypath, models[model_version])
time.sleep(trace_id) # make sure this doesn't crash the script if multiple instances of the model are submitted simultaneously
if not os.path.exists(thispath):
    try:
        os.mkdir(thispath)
    except:
        pass

# and... go!
starttime = time.time()
run_model(mypath, models[model_version], trace_id)
elapsed = time.time() - starttime
print( "Elapsed time: %f seconds\n" %elapsed )

# and plot
# import HDDM_plotOutput
# HDDM_plotOutput.plot_model(mypath, models[model_version], trace_id)
