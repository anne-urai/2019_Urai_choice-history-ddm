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
        default = 1,
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

def run_model(mypath, model_name, trace_id, nr_samples=50000):

    import os
    import hddm

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename
    modelExists     = os.path.isfile(model_filename)

    if not modelExists:
        # get the csv
        mydata = hddm.load_csv(os.path.join(mypath, '2ifc_data_hddm.csv'))

        # specify the model
        if model_name == 'stimcoding':
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr']},
                p_outlier=0.05)

        elif model_name == 'prevresp_z':
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'z':['prevresp']},
                p_outlier=0.05)

        elif model_name == 'prevresp_dc':
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'dc':['prevresp']},
                p_outlier=0.05)

        elif model_name == 'prevresp_prevrt_z':
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'z':['prevresp', 'prevrt']},
                p_outlier=0.05)

        elif model_name == 'prevresp_prevrt_dc':
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'dc':['prevresp', 'prevrt']},
                p_outlier=0.05)

        elif model_name == 'prevresp_prevpupil_z':
            # take only rows without NaN
            m = hddm.HDDMStimCoding(mydata.dropna(), stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'z':['prevresp', 'prevpupil']},
                p_outlier=0.05)

        elif model_name == 'prevresp_prevpupil_dc':
            # take only rows without NaN
            m = hddm.HDDMStimCoding(mydata.dropna(), stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v':['sessionnr'], 'dc':['prevresp', 'prevpupil']},
                p_outlier=0.05)

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

# ============================================ #
# run one model per job
# ============================================ #

# find path depending on location
import os, time
usr = os.environ.get('USER')

if usr in ['anne']: #local
    mypath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/Data/HDDM'
if usr in ['aurai']: #uke cluster
    mypath = '/home/aurai/Data/MEG-PL/Data/HDDM'
if usr in ['aeurai']: #cartesius
    mypath = '/home/aeurai/neurodec/Data/MEG-PL/HDDM'
    # to avoid errors when plotting on cartesius
    # http://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
    import matplotlib
    matplotlib.use('Agg') # to still plot even when no display is defined

# which model are we running at the moment?
models = {0: 'stimcoding',
    1: 'prevresp_z',
    2: 'prevresp_dc',
    3: 'prevresp_prevrt_z',
    4: 'prevresp_prevrt_dc',
    5: 'prevresp_prevpupil_z',
    6: 'prevresp_prevpupil_dc'}

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
import HDDM_plotOutput
HDDM_plotOutput.plot_model(mypath, models[model_version], trace_id)
