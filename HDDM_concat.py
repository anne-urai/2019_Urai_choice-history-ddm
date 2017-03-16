#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2016
takes input arguments from stopos
Important: on Cartesius, call module load pytho/2.7.9 before running
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
        default = range(1,2),
        type = "int",
        help = "Which dataset, see below" )
parser.add_option ( "-v", "--version",
        default = range(0,5),
        type = "int",
        help = "Version of the model to run" )

opts,args       = parser.parse_args()
model_version   = opts.version
d               = opts.dataset

# ============================================ #
# define the function that will do the work
# ============================================ #

def concat_models(mypath, model_name):

    import os, hddm, time, kabuki

    # ============================================ #
    # APPEND CHAINS
    # ============================================ #

    allmodels = []
    print "appending models"
    for trace_id in range(15): # 15 models were run
        model_filename              = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
        print(model_filename)
        modelExists                 = os.path.isfile(model_filename)
        assert modelExists == True # if not, this model has to be rerun

        starttime = time.time()
        thism                       = hddm.load(model_filename)

        # now append
        allmodels.append(thism)
        elapsed = time.time() - starttime
        print( "Elapsed time: %f seconds." %elapsed )

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
    m = kabuki.utils.concat_models(allmodels)
    m.save(os.path.join(mypath, model_name, 'modelfit-combined.model')) # save the model to disk

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

    # can then plot full posteriors and compute p-values in Matlab

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

datasets = {0: 'RT_RDK', 1: 'MEG-PL'}

# recode
if isinstance(d, int):
    d = range(d,d+1) # makes a list out of an integer
if isinstance(model_version, int):
    model_version = range(model_version, model_version+1)

for dx in d:
    # find path depending on location and dataset
    import os, time
    mypath = os.path.realpath(os.path.expanduser('~/neurodec/Data/%s/HDDM'%datasets[dx]))

    for vx in model_version:
        # and... go!
        starttime = time.time()
        concat_models(mypath, models[vx])
        elapsed = time.time() - starttime
        print( "Elapsed time: %f seconds\n" %elapsed )

    # and plot
    # import HDDM_plotOutput
    # HDDM_plotOutput.plot_model(mypath, models[model_version], trace_id)
