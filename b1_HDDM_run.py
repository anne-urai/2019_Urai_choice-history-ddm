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
import os, fnmatch
import corner
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
        default = range(0,7),
        type = "int",
        help = "Which dataset, see below" )
parser.add_option ( "-v", "--version",
        default = range(0,11),
        type = "int",
        help = "Version of the model to run" )
parser.add_option ( "-i", "--trace_id",
        default = 14,
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
    if not (os.path.isfile(os.path.join(mypath, model_name, 'modelfit-md14.model'))) and  (os.path.isfile(os.path.join(mypath, model_name, 'modelfit-combined.model'))):
        print os.path.join(mypath, model_name, 'modelfit-combined.model')
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

        if len(allmodels) == 0:
            return allmodels

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
        # SAVE THE FULL MODEL
        # ============================================ #

        print "concatenated models"
        m.save(os.path.join(mypath, model_name, 'modelfit-combined.model')) # save the model to disk

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

def cornerplot(mypath, datasetname, modelname):

    # ============================================ #
    # corner plot for parameter recovery
    # ============================================ #

    m = hddm.load(os.path.join(mypath, modelname, 'modelfit-combined.model'))
    print os.path.join(mypath, modelname, 'modelfit-combined.model')

    # dictionary with all parameters that were fit
    params_of_interest_0 = list(m.get_group_nodes().index)
    # dont plot the std
    params_of_interest_0 = [x for x in params_of_interest_0 if not 'std' in x]
    params_of_interest_0 = [x for x in params_of_interest_0 if not 'sv' in x]

    traces_0 = []
    for p in range(len(params_of_interest_0)):
        traces_0.append(m.nodes_db.node[params_of_interest_0[p]].trace.gettrace())

    fig = corner.corner(np.array(traces_0).T, color='b', labels=params_of_interest_0, show_titles=True, **{'lw':1})
    try:
        fig.savefig(os.path.join('/nfs/aeurai/HDDM/summary/figures',  'corner_%s_%s.pdf' %(datasetname,modelname)))
    except:
        print('cannot save figure')

    df0 = pd.DataFrame(np.array(traces_0).T[:,:len(params_of_interest_0)], columns=params_of_interest_0)
    fig = corner.corner(df0, color='b', **{'lw':1})

    # now add regression lines and stats, from JW
    for i, j in zip(*np.triu_indices_from(np.zeros((len(params_of_interest_0),len(params_of_interest_0))), 1)):
        # add titles:
        r0, p0 = sp.stats.pearsonr(df0.iloc[:,i], df0.iloc[:,j])
        fig.axes[(j*len(params_of_interest_0))+i].set_title('r={}; p={}'.format(round(r0, 3), round(p0, 3),))
        # add regression lines:
        x_line = np.linspace(fig.axes[(j*len(params_of_interest_0))+i].axis()[0], fig.axes[(j*len(params_of_interest_0))+i].axis()[1], 100)
        (m,b) = sp.polyfit(df0.iloc[:,i], df0.iloc[:,j],1)
        regression_line = sp.polyval([m,b],x_line)
        if p0 < 0.05:
            fig.axes[(j*len(params_of_interest_0))+i].plot(x_line, regression_line, color='r', zorder=3)
        else:
            fig.axes[(j*len(params_of_interest_0))+i].plot(x_line, regression_line, color='b', zorder=3)

    # sns.despine(offset=0, trim=True)
    plt.tight_layout()
    try:
        fig.savefig(os.path.join('/nfs/aeurai/HDDM/summary/figures',  'corner_%s_%s.pdf' %(datasetname,modelname)))
    except:
        print('cannot save figure')

# ============================================ #
# PREPARE THE ACTUAL MODEL FITS
# ============================================ #

# which model are we running at the moment?
models = ['stimcoding_nohist', # 0
    'stimcoding_dc_prevresp', # 1
    'stimcoding_z_prevresp', # 2
    'stimcoding_dc_z_prevresp', # 3
    'stimcoding_dc_prevcorrect', # 4
    'stimcoding_z_prevcorrect', # 5
    'stimcoding_dc_z_prevcorrect', # 6
    'regress_dc_z_prevresp', # 7
    'regress_dc_z_prevresp_prevrt', # 8
    'regress_dc_z_prev2resp', # 9
    'regress_dc_z_prev3resp', # 10
    'regress_dc_z_prevresp_prevstim_prevrt_prevpupil', # 11
    'stimcoding_dc_z_prevresp_pharma', #12
    'stimcoding_dc_z_prevresp_sessions', # 13
    'regress_dc_z_visualgamma', #14
    'regress_dc_z_motorslope', #15
    'regress_dc_z_motorstart', #16
    'stimcoding_sz_nohist', # 17
    'stimcoding_sz_dc_prevresp', # 18
    'stimcoding_sz_z_prevresp', # 19
    'stimcoding_sz_dc_z_prevresp', # 20
    ]

datasets = ['RT_RDK', # 0
    'MEG', # 1
    'NatComm', # 2
    'Anke_merged', # 3
    'JW_yesno', # 4
    'Bharath_fMRI', # 5
    'Murphy', # 6
    'MEG_MEGdata'] # 7
    
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
    elif 'anne' in usr:
        mypath = os.path.realpath(os.path.expanduser('~/Data/HDDM/%s'%datasets[dx]))

    for vx in model_version:
        time.sleep(trace_id) # to avoid different jobs trying to make the same folder

        # make a folder for the outputs, combine name and time
        thispath = os.path.join(mypath, models[vx])
        if not os.path.exists(thispath):
            os.mkdir(thispath)

        if runMe == 1:

            starttime = time.time()
            model_filename = os.path.join(mypath, models[vx], 'modelfit-md%d.model'%trace_id)

            # get the csv file for this dataset
            filename    = fnmatch.filter(os.listdir(mypath), '*.csv')
            mydata      = hddm.load_csv(os.path.join(mypath, filename[0]))

            # correct a weirdness in Anke's data
            if 'transitionprob' in mydata.columns:
                mydata.transitionprob = mydata.transitionprob * 100;
                mydata.transitionprob = mydata.transitionprob.round();

            # get the model specification, pass data
            m = make_model(mypath, mydata, models[vx], trace_id)

            # now sample and save
            if os.path.exists(model_filename):
                pass # skip if this model i has been run
            elif os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')) and not os.path.exists(model_filename):
                pass # skip if this model has been concatenated
            else:
                # only run if this hasnt been done, and there is no concatenated master model present
                run_model(m, mypath, models[vx], trace_id, n_samples)
            elapsed = time.time() - starttime
            print( "Elapsed time for %s, %s, %d samples: %f seconds\n" %(models[vx], datasets[dx], n_samples, elapsed))

            # ================================================= #
            # important, concat after running to save disk space
            # ================================================= #

            if trace_id == 14 and not os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')):
                # https://stackoverflow.com/questions/35795452/checking-if-a-list-of-files-exists-before-proceeding
                filelist = []
                for t in range(15):
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

            # make corner plot
            if trace_id == 14 and os.path.exists(os.path.join(mypath, models[vx], 'modelfit-combined.model')):
                cornerplot(mypath, datasets[dx], models[vx])

        elif runMe == 2:

            # ============================================ #
            # POSTERIOR PREDICTIVES TO ASSESS MODEL FIT
            # ============================================ #

            starttime = time.time()
            print "computing ppc"

            # specify how many samples are needed
            m = hddm.load(os.path.join(mypath, models[vx], 'modelfit-combined.model'))
            print os.path.join(mypath, models[vx], 'modelfit-combined.model')
            ppc = hddm.utils.post_pred_gen(m, append_data=True, samples=100)

            # make the csv smaller, save disk space
            ppc = ppc[['rt','rt_sampled', 'stimulus', 'response']]

            # save as pandas dataframe
            ppc.to_csv(os.path.join(mypath, models[vx], 'ppq_data.csv'), index=True)
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
            
            subj_params = []
            bic         = []
            for subj_idx, subj_data in mydata.groupby('subj_idx'):
                m_subj = make_model(mypath, subj_data, models[vx], trace_id)
                thismodel = m_subj.optimize('gsquare')
                thismodel.update({'subj_idx':subj_idx}) # keep original subject number
                subj_params.append(thismodel)
                bic.append(m_subj.bic_info)
                
            params = pd.DataFrame(subj_params)
            params.to_csv(os.path.join(mypath, models[vx], 'chisquare.csv'))
            bic = pd.DataFrame(bic)
            bic.to_csv(os.path.join(mypath, models[vx], 'BIC.csv'))
           
