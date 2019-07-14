#!/usr/bin/env python
# encoding: utf-8

"""
% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

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
import os, fnmatch
import pandas as pd
import scipy as sp

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
        default = range(26),
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 29,
        type = "int",
        help = "Which trace to run, usually 0-60" )
parser.add_option ( "-s", "--samples",
        default = 50,
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
    m.sample(n_samples, burn=n_samples/2, thin=3, db='pickle',
        dbname=os.path.join(mypath, model_name, 'modelfit-md%d.db'%trace_id))
    m.save(os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)) # save the model to disk

    # ============================================ #
    # save the output values
    # ============================================ #

    # save the DIC for this model
    text_file = open(os.path.join(mypath, model_name, 'DIC-md%d.txt'%trace_id), 'w')
    text_file.write("Model {}: {}\n".format(trace_id, m.dic))
    text_file.close()

    # save the other model comparison indices
    df = dict()
    df['dic_original'] = [m.dic]
    df['aic'] = [aic(m)]
    df['bic'] = [bic(m)]
    df2 = pd.DataFrame(df)
    df2.to_csv(os.path.join(mypath, models[vx], 'model_comparison_md%d.csv'%trace_id))


def concat_models(mypath, model_name):
    
    nchains = 30

    # CHECK IF COMBINED MODEL EXISTS
    if os.path.isfile(os.path.join(mypath, model_name, 'modelfit-combined.model')):
        print os.path.join(mypath, model_name, 'modelfit-combined.model')

    else:
        # ============================================ #
        # APPEND MODELS
        # ============================================ #

        allmodels = []
        print ("appending models for %s" %model_name)
        for trace_id in range(nchains): # how many chains were run?
            model_filename        = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
            modelExists           = os.path.isfile(model_filename)

            if modelExists == True: # if not, this model has to be rerun
                print model_filename
                thism                 = hddm.load(model_filename)
                allmodels.append(thism) # now append into a list

        # ============================================ #
        # CHECK CONVERGENCE
        # ============================================ #

        if len(allmodels) == 0:
            return allmodels

        try:
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
        except:
            pass

        # now actually concatenate them, see email Gilles
        m = kabuki.utils.concat_models(allmodels)

        # ============================================ #
        # SAVE THE FULL MODEL
        # ============================================ #

        print "concatenated models"
        m.save(os.path.join(mypath, model_name, 'modelfit-combined.model')) # save the model to disk

        # DELETE FILES to save space
        print "deleting separate chains..."
        for fl in glob.glob(os.path.join(mypath, models[vx], 'modelfit-md*.model')):
                print(fl)
                os.remove(fl)
        for fl in glob.glob(os.path.join(mypath, models[vx], 'modelfit-md*.db')):
            if not '-md0.db' in fl:
                print(fl)
                os.remove(fl)

        # ============================================ #
        # SAVE POINT ESTIMATES
        # ============================================ #

        print "saving stats"
        results = m.gen_stats() # point estimate for each parameter and subject
        results.to_csv(os.path.join(mypath, model_name, 'results-combined.csv'))

        # save the DIC for this model
        text_file = open(os.path.join(mypath, model_name, 'DIC-combined.txt'), 'w')
        text_file.write("Combined model: {}\n".format(m.dic))
        text_file.close()

        # ============================================ #
        # SAVE TRACES
        # ============================================ #

        print "saving traces"
        # get the names for all nodes that are available here
        group_traces = m.get_group_traces()
        group_traces.to_csv(os.path.join(mypath, model_name, 'group_traces.csv'))

        all_traces = m.get_traces()
        all_traces.to_csv(os.path.join(mypath, model_name, 'all_traces.csv'))


        # ============================================ #
        # CONCATENATE MODEL COMPARISON
        # ============================================ #

        # average model comparison values across chains
        print('concatenating model comparison')

        nchains = 30
        for trace_id in range(nchains): # how many chains were run?
            filename = os.path.join(mypath, models[vx], 'model_comparison_md%d.csv'%trace_id)
            df = pd.read_csv(filename)

            if trace_id == 0:
                df2 = df
            else:
                df2 = df2.append(df, ignore_index=True)

        # average over chains
        df3 = df2.mean()
        df3 = df2.describe().loc[['mean']]
        df3.to_csv(os.path.join(mypath, models[vx], 'model_comparison.csv'))

        for fl in glob.glob(os.path.join(mypath, models[vx], 'model_comparison_md*.csv')):
            print(fl)
            os.remove(fl)
            

# ============================================ #
# also compute BIC, AIC
# from https://groups.google.com/forum/#!searchin/hddm-users/bic%7Csort:date/hddm-users/Bo2vUcpR008/RLRpL0faptAJ
# ============================================ #

def aic(self):
    k = len(self.get_stochastics())
    logp = sum([x.logp for x in self.get_observeds()['node']])  
    return 2 * k - 2 * logp

def bic(self):
    k = len(self.get_stochastics())
    n = len(self.data)
    logp = sum([x.logp for x in self.get_observeds()['node']])
    return -2 * logp + k * np.log(n)

# ============================================ #
# PREPARE THE ACTUAL MODEL FITS
# ============================================ #

# SPECIFY THE MODELS TO RUN HERE - SEE hddm_models.py for specifications
models = ['regress_nohist', #0
'regress_dc_lag1',
'regress_z_lag1',
'regress_dcz_lag1',
'regress_dc_lag2',
'regress_z_lag2',
'regress_dcz_lag2',
'regress_dc_lag3',
'regress_z_lag3',
'regress_dcz_lag3',
'regress_dc_lag4', #14
'regress_z_lag4',
'regress_dcz_lag4',
'regress_dc_lag5',
'regress_z_lag5',
'regress_dcz_lag5',
'regress_dc_lag6', #20
'regress_z_lag6',
'regress_dcz_lag6',
'regress_dc_lag7',
'regress_z_lag7',
'regress_dcz_lag7', # 2
'stimcoding_nohist', # 0
'stimcoding_dc_prevresp', #1
'stimcoding_z_prevresp', #2
'stimcoding_dc_z_prevresp', #3
'stimcoding_dc_z_prevresp_st', #4
'stimcoding_dc_z_prevresp_pharma', #5
'stimcoding_dc_z_prevcorrect',#6
'stimcoding_prevcorrect',#7
'stimcoding_dc_z_prev2resp', #
'stimcoding_dc_prevcorrect', #11
'regress_dc_z_visualgamma',  #12
'regress_dc_z_motorstart', #13
'regress_dc_z_prevresp_visualgamma', #14
'regress_dc_z_prevresp_motorstart', #15
'stimcoding_nohist_stcoh',  # `10
'stimcoding_dc_prevresp_stcoh',  # 
'stimcoding_z_prevresp_stcoh',  # 
'stimcoding_dc_z_prevresp_stcoh']  # ] # 50

models = ['stimcoding_nohist',
          'stimcoding_dc_prevresp',
          'stimcoding_z_prevresp',
          'stimcoding_dc_z_prevresp',
          'stimcoding_dc_prevstim',
          'stimcoding_z_prevstim',
          'stimcoding_dc_z_prevresp']

datasets = ['Murphy', 'JW_yesno', 'JW_PNAS', 'NatComm', 'MEG', 
    'Anke_MEG_neutral', 'Anke_MEG_transition', 'Anke_MEG_transition_no81', 
    'MEG_MEGdata',  'NatComm_rescaled']

datasets = ['Murphy', 'JW_yesno', 'JW_PNAS', 'NatComm', 'MEG', 'Anke_MEG_transition']
datasets = ['Anke_MEG_blocks']

# recode
if isinstance(d, int):
    d = range(d,d+1) # makes a list out of an integer
if isinstance(model_version, int):
    model_version = range(model_version, model_version+1)

for dx in d:

    # find path depending on location and dataset
    usr = os.environ['USER']
    if 'aeurai' in usr:
        mypath = os.path.realpath(os.path.expanduser('/nfs/aeurai/HDDM/%s'%datasets[dx]))
        # LISA PROJECT SPACE ENDED, USE HOME SPACE
        mypath = os.path.realpath(os.path.expanduser('/home/aeurai/Data/HDDM/%s'%datasets[dx]))
    elif 'anne' in usr:
        mypath = os.path.realpath(os.path.expanduser('~/Data/HDDM/%s'%datasets[dx]))

    for vx in model_version:
        time.sleep(trace_id) # to avoid different jobs trying to make the same folder

        # make a folder for the outputs, combine name and time
        thispath = os.path.join(mypath, models[vx])
        if not os.path.exists(thispath):
            os.mkdir(thispath)
    
        if runMe == 1:

            # get the csv file for this dataset
            filename    = fnmatch.filter(os.listdir(mypath), '*.csv')
            mydata      = hddm.load_csv(os.path.join(mypath, filename[0]))

            # round up
            if 'transitionprob' in mydata.columns:
                mydata.transitionprob = mydata.transitionprob * 100;
                mydata.transitionprob = mydata.transitionprob.round();

            starttime = time.time()
            model_filename = os.path.join(mypath, models[vx], 'modelfit-md%d.model'%trace_id)

            # ============================================ #
            # DECIDE WHAT TO DO
            # ============================================ #

            # now sample and save
            if not os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')) and os.path.exists(model_filename):
                pass # this model has been run but the job isn't finished
            elif os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')) and os.path.exists(model_filename):
                
                # there is a concatenated model, but this file still remains - delete!
                # DELETE FILES to save space
                print "deleting separate chains..."
                for fl in glob.glob(os.path.join(mypath, models[vx], 'modelfit-md*.model')):
                        print(fl)
                        os.remove(fl)
                for fl in glob.glob(os.path.join(mypath, models[vx], 'modelfit-md*.db')):
                    if not '-md0.db' in fl:
                        print(fl)
                        os.remove(fl)

            elif os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')) and not os.path.exists(model_filename):
                print('model already run and concatenated, exiting...')
                pass # model has been concatenated, cleanup done

            elif models[vx] == 'stimcoding_dc_z_prevresp_pharma' and not 'drug' in mydata.columns:
                continue # makes no sense, don't run
            elif 'stcoh' in models[vx] and not 'coherence' in mydata.columns:
                continue # makes no sense, don't run

            else:
                print( "Starting to run %s, %s, %d samples" %(models[vx], datasets[dx], n_samples))
                # get the model specification, pass data
                m = make_model(mypath, mydata, models[vx], trace_id)
                # only run if this hasnt been done, and there is no concatenated master model present
                run_model(m, mypath, models[vx], trace_id, n_samples)
                
            elapsed = time.time() - starttime
            print( "Elapsed time for %s, %s, %d samples: %f seconds\n" %(models[vx], datasets[dx], n_samples, elapsed))

            # ================================================= #
            # important, concat after running to save disk space
            # ================================================= #

            if trace_id == 29 and not os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')): # and not os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')):
                # https://stackoverflow.com/questions/35795452/checking-if-a-list-of-files-exists-before-proceeding
                filelist = []
                for t in range(30):
                    filelist.append(os.path.join(mypath, models[vx], 'modelfit-md%d.model'%t))

                print filelist
                # wait until all the files are present
                while True:
                    if all([os.path.isfile(f) for f in filelist]):
                        break
                    else: # wait
                        print("waiting for files")
                        # raise ValueError('Not all files present')
                        time.sleep(60)

                # concatenate the different chains, will save disk space
                concat_models(mypath, models[vx])

        elif runMe == 2:

            # ============================================ #
            # POSTERIOR PREDICTIVES TO ASSESS MODEL FIT
            # ============================================ #

            starttime = time.time()
            print "computing ppc"

            # specify how many samples are needed
            m = hddm.load(os.path.join(mypath, models[vx], 'modelfit-combined.model'))
            print os.path.join(mypath, models[vx], 'modelfit-combined.model')
            if 'MEG' in datasets[dx]:
                nsmp = 50
            else:
                nsmp = 100

            ppc = hddm.utils.post_pred_gen(m, append_data=True, samples=nsmp)

            # make the csv smaller, save disk space
            savecols = list(set(ppc.columns) & set(['rt','rt_sampled', 'response_sampled',
                        'index', 'stimulus', 'response', 'prevresp', 'subj_idx',
                        'transitionprob', 'coherence', 'prevcorrect']))
            ppc = ppc[savecols]

            # save as pandas dataframe
            ppc.to_csv(os.path.join(mypath, models[vx], 'ppc_data.csv'), index=True)
            elapsed = time.time() - starttime
            print( "Elapsed time for %s %s, PPC: %f seconds\n" %(models[vx], datasets[dx], elapsed))

        elif runMe == 3:

            # ============================================ #
            # QUANTILE OPTIMISATION
            # http://ski.clps.brown.edu/hddm_docs/howto.html#run-quantile-opimization
            # ============================================ #

            # get the csv file for this dataset
            filename    = fnmatch.filter(os.listdir(mypath), '*.csv')
            mydata      = hddm.load_csv(os.path.join(mypath, filename[0]))
            # mydata      = mydata[mydata.rt > 0.25] # remove superfast responses

            subj_params = []
            bic         = []
            for subj_idx, subj_data in mydata.groupby('subj_idx'):
                m_subj    = make_model(mypath, subj_data, models[vx], trace_id)
                # m_subj.find_starting_values() # this may help the fits
                # thismodel = m_subj.optimize('gsquare')
                thismodel = m_subj.optimize('gsquare', quantiles=[0.1, 0.3, 0.5, 0.7, 0.9], n_runs=5)
                thismodel.update({'subj_idx':subj_idx}) # keep original subject number
                subj_params.append(thismodel)
                bic.append(m_subj.bic_info)

            params = pd.DataFrame(subj_params)
            params.to_csv(os.path.join(mypath, models[vx], 'Gsquare.csv'))
            bic = pd.DataFrame(bic)
            bic.to_csv(os.path.join(mypath, models[vx], 'BIC.csv'))







