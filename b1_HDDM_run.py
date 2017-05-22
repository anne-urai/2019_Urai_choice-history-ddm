#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2017
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

# to avoid errors when plotting on cartesius
# http://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
import matplotlib
matplotlib.use('Agg') # to still plot even when no display is defined
import matplotlib.pyplot as plt
from IPython import embed as shell
import numpy as np

import warnings
warnings.filterwarnings('ignore')

# get the model specification here
from hddm_models import make_model
import os, hddm, time, kabuki, glob
from math import ceil

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
        default = range(0,6),
        type = "int",
        help = "Which dataset, see below" )
parser.add_option ( "-v", "--version",
        default = range(0,25),
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 0,
        type = "int",
        help = "Which trace to run, usually 0-60" )
parser.add_option ( "-s", "--samples",
        default = 100,
        type = "int",
        help = "How many samples to use" )

opts,args       = parser.parse_args()
model_version   = opts.version
d               = opts.dataset
trace_id        = opts.trace_id
runMe           = opts.run
n_samples       = opts.samples

def run_model(m, mypath, model_name, trace_id, n_samples):

    # ============================================ #
    # do the actual sampling
    # ============================================ #

    print "finding starting values"
    m.find_starting_values() # this should help the sampling

    print "begin sampling"
    m.sample(n_samples, burn=n_samples/2, thin=2, db='pickle',
        dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
    # specify a certain backend? pickle?
    m.save(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)) # save the model to disk

    # ============================================ #
    # save the output values
    # ============================================ #

    # save the DIC for this model
    text_file = open(os.path.join(mypath, model_name, 'DIC-md%d.txt'%trace_id), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

def concat_models(mypath, model_name):

    # CHECK IF COMBINED MODEL EXISTS
    if not (os.path.isfile(os.path.join(mypath, model_name, 'modelfit-md14.model'))) and (os.path.isfile(os.path.join(mypath, model_name, 'modelfit-combined.model'))):
        m = hddm.load(os.path.join(mypath, model_name, 'modelfit-combined.model'))
    else:

        # ============================================ #
        # APPEND MODELS
        # ============================================ #

        allmodels = []
        print ("appending models for %s" %model_name)
        for trace_id in range(15): # how many chains were run?
            model_filename        = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
            modelExists           = os.path.isfile(model_filename)
            if modelExists == True: # if not, this model has to be rerun
                print model_filename
                thism                 = hddm.load(model_filename)
                allmodels.append(thism) # now append into a list

        # ============================================ #
        # CHECK CONVERGENCE
        # ============================================ #

        gr = hddm.analyze.gelman_rubin(allmodels)

        # save
        text_file = open(os.path.join(mypath, model_name, 'gelman_rubin.txt'), 'w')
        for p in gr.items():
            text_file.write("%s,%s\n" % p)
            # print a warning when non-convergence is detected
            # Values should be close to 1 and not larger than 1.02 which would indicate convergence problems.
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3731670/
            if abs(p[1]-1) > 0.02:
                print "non-convergence found, %s:%s" %p
        text_file.close()
        print "written gelman rubin stats to file"

        # now actually concatenate them, see email Gilles
        m = kabuki.utils.concat_models(allmodels)
        print "concatenated models"
        m.save(os.path.join(mypath, model_name, 'modelfit-combined.model')) # save the model to disk
        # m.save(os.path.join(mypath, model_name, 'modelfit-combined.pickle')) # save the model to disk

        # ============================================ #
        # DELETE FILES to save space
        # ============================================ #

        if len(allmodels) == 15:
            print "deleting separate chains"
            for fl in glob.glob(os.path.join(mypath, model_name, 'modelfit-md*.model')):
                    os.remove(fl)
            for fl in glob.glob(os.path.join(mypath, model_name, 'modelfit-md*.db')):
                if not '-md0.db' in fl:
                    os.remove(fl)
        else:
            print "not deleting individual model chains"

    # ============================================ #
    # POSTERIOR PREDICTIVE PLOTS
    # ============================================ #

    if datasets[dx] == 'MEG':
        num_subj = 65
    if datasets[dx] == 'RT_RDK':
        num_subj = 25
    if datasets[dx] in ['Anke_serial', 'Anke_neutral', 'Anke_repetitive', 'Anke_alternating']:
        num_subj = 27
    if datasets[dx] == 'NatComm':
        num_subj = 27

    figpath = os.path.join(mypath, model_name, 'figures-concat')
    if not os.path.exists(figpath):
        os.mkdir(figpath)

    m.plot_posterior_predictive(save=True, path=figpath, format='pdf',
        columns=ceil(num_subj/6))
    plt.close('all') # to avoid warnings
    print "plotted posterior predictive RT distributions"

    # plot the traces and posteriors for each parameter
    m.plot_posteriors(save=True, path=figpath, format='pdf')
    plt.close('all') # to avoid warnings
    print "plotted traces and autocorrelation"

    # ======================================================================================== #
    # plot my own version of the posterior predictive with overlaid conditions
    # ======================================================================================== #

    if 'stimcoding' in model_name:
        print "exporting posterior predictives"
        from kabuki.analyze import _plot_posterior_pdf_node

        if not os.path.exists(os.path.join(mypath, model_name, 'preds')):
            os.mkdir(os.path.join(mypath, model_name, 'preds'))
        observeds = m.get_observeds()

        # Plot different conditions (new figure for each)
        for tag, nodes in observeds.groupby('tag'):
            # retrieve individual subjects
            for subj_i, (node_name, bottom_node) in enumerate(nodes.iterrows()):
                fig = plt.figure()
                ax  = fig.add_subplot(2,2,1)
                #    fig.suptitle(utils.pretty_tag(tag), fontsize=8)
                fig.subplots_adjust(top=0.9, hspace=.4, wspace=.3)

                if not hasattr(bottom_node['node'], 'pdf'):
                    continue # skip nodes that do not define the required_method
                y = _plot_posterior_pdf_node(bottom_node['node'], ax,
                    value_range=np.linspace(-3,3,100))
                rtvals = bottom_node['node'].value.values

                # save this figure
                fig.savefig(os.path.join(mypath, model_name, 'preds', 'ppq_%s_subj%d.pdf'%(str(tag),subj_i)))
                plt.close()

                # now save to a file so that I can plot it in matlab
                np.savetxt(os.path.join(mypath, model_name, 'preds', 'ppq_y_%s_subj%d.csv'%(str(tag),subj_i)), y, delimiter=",")
                np.savetxt(os.path.join(mypath, model_name, 'preds', 'ppq_rt_%s_subj%d.csv'%(str(tag),subj_i)), rtvals, delimiter=",")

    # ============================================ #
    # SAVE POINT ESTIMATES
    # ============================================ #

    results = m.gen_stats() # point estimate for each parameter and subject
    results.to_csv(os.path.join(mypath, model_name, 'results-combined.csv'))

    # save the DIC for this model
    text_file = open(os.path.join(mypath, model_name, 'DIC-combined.txt'), 'w')
    text_file.write("Combined model: {}\n".format(m.dic))
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
models = ['stimcoding_nohist', # 0
    'regress_nohist', # 1
    'stimcoding_dc_prevresp', # 2
    'stimcoding_z_prevresp', # 3
    'stimcoding_dc_z_prevresp', # 4
    'stimcoding_dc_prevresp_prevstim', # 5
    'stimcoding_z_prevresp_prevstim', # 6
    'stimcoding_dc_z_prevresp_prevstim', # 7
    'regress_dc_prevresp', # 8
    'regress_z_prevresp', # 9
    'regress_dc_z_prevresp', # 10
    'regress_dc_prevresp_prevstim', # 11
    'regress_z_prevresp_prevstim', # 12
    'regress_dc_z_prevresp_prevstim', # 13
    'regress_dc_prevresp_prevstim_vasessions', # 14
    'regress_dc_prevresp_prevstim_prevpupil', # 15
    'regress_dc_prevresp_prevstim_prevrt', # 16
    'regress_dc_prevresp_prevstim_prevrt_prevpupil', # 17
    'regress_dc_z_prevresp_prevstim_vasessions', # 18
    'regress_dc_z_prevresp_prevstim_prevpupil', # 18
    'regress_dc_z_prevresp_prevstim_prevrt', # 20
    'regress_dc_z_prevresp_prevstim_prevrt_prevpupil', # 21
    'regress_dc_prevresp_prevstim_vasessions_prevrespsessions', # 22
    'regress_dc_prevresp_prevstim_vasessions_prevpupil', # 23
    'regress_dc_prevresp_prevstim_vasessions_prevrt', # 24
    'regress_dc_prevresp_prevstim_vasessions_prevrt_prevpupil', # 25
    'regress_dc_prev2resp_prev2stim', # 26
    'regress_dc_prev3resp_prev3stim'] # 27

datasets = ['RT_RDK', # 0
    'MEG', # 1
    'NatComm', # 2
    'Anke_neutral', # 3
    'Anke_repetitive', # 4
    'Anke_alternating', # 5
    'Anke_serial'] # 6

# recode
if isinstance(d, int):
    d = range(d,d+1) # makes a list out of an integer
if isinstance(model_version, int):
    model_version = range(model_version, model_version+1)

for dx in d:

    # find path depending on location and dataset
    mypath = os.path.realpath(os.path.expanduser('~/Data/%s/HDDM'%datasets[dx]))

    for vx in model_version:
        time.sleep(trace_id) # to avoid different jobs trying to make the same folder

        # make a folder for the outputs, combine name and time
        thispath = os.path.join(mypath, models[vx])
        if not os.path.exists(thispath):
            os.mkdir(thispath)

        if runMe == True:
            starttime = time.time()
            model_filename = os.path.join(mypath, models[vx], 'modelfit-md%d.model'%trace_id)
            # get the model specification
            m = make_model(mypath, models[vx], trace_id)

            # now sample and save
            if os.path.exists(model_filename):
                pass
            elif os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')) and not os.path.exists(model_filename):
                pass
            else:
                # only run if this hasnt been done, and there is no concatenated master model present
                run_model(m, mypath, models[vx], trace_id, n_samples)
            elapsed = time.time() - starttime
            print( "Elapsed time for %s, %s, %d samples: %f seconds\n" %(models[vx], datasets[dx], n_samples, elapsed))

        else: # concatenate the different chains
            concat_models(mypath, models[vx])
