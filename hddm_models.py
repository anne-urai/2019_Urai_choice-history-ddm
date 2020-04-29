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
import numpy as np
import os, fnmatch
import hddm
from IPython import embed as shell
import pandas as pd

# prepare link function for the regression models
def z_link_func(x):
    return 1 / (1 + np.exp(-(x.values.ravel())))

def balance_designmatrix(mydata):
    # remove subjects who did not do all conditions
    for i, sj in enumerate(mydata.subj_idx.unique()):
        sessions = mydata[mydata.subj_idx == sj].session.unique()
        if len(sessions) < len(mydata.session.unique()):
            mydata = mydata[mydata.subj_idx != sj] # drop this subject
    return mydata

def recode_4stimcoding(mydata):
    # split into coherence and stimulus identity
    mydata['coherence'] = mydata.stimulus.abs()
    mydata.stimulus     = np.sign(mydata.stimulus)
    # for stimcoding, the two identities should be 0 and 1
    mydata.ix[mydata['stimulus']==-1,'stimulus'] = 0
    if len(mydata.stimulus.unique()) != 2:
        raise ValueError('Stimcoding needs 2 stimulus types')

    # also create a binary prevcorrect
    mydata['prevcorrect']     = mydata.prevresp
    mydata.prevcorrect[mydata.prevresp != mydata.prevstim] = 0
    mydata.prevcorrect[mydata.prevresp == mydata.prevstim] = 1

    try:
        # also create a binary prevcorrect
        mydata['prev2correct']     = mydata.prevresp
        mydata.prev2correct[mydata.prev2resp != mydata.prev2stim] = 0
        mydata.prev2correct[mydata.prev2resp == mydata.prev2stim] = 1
        # also create a binary prevcorrect
        mydata['prev3correct']     = mydata.prevresp
        mydata.prev3correct[mydata.prev3resp != mydata.prev3stim] = 0
        mydata.prev3correct[mydata.prev3resp == mydata.prev3stim] = 1
    except:
        pass

    return mydata

def add_morelags(mydata):

    # make sure stimulus and response are [-1,1]
    response = np.sign(mydata.response - 0.1)
    stimulus = np.sign(mydata.stimulus)

    for lag in range(16):
        if not 'prev%dresp'%lag in mydata.columns:
            mydata['prev%dresp'%lag] = np.roll(response, lag)

        if not 'prev%dstim'%lag in mydata.columns:
            mydata['prev%dstim'%lag] = np.roll(stimulus, lag)

        # recode into previous correct and error
        mydata['prev%dresp_correct'%lag] = mydata['prev%dresp'%lag]
        mydata['prev%dresp_correct'%lag][mydata['prev%dresp'%lag] != mydata['prev%dstim'%lag]] = 0
        mydata['prev%dresp_error'%lag] = mydata['prev%dresp'%lag]
        mydata['prev%dresp_error'%lag][mydata['prev%dresp'%lag] == mydata['prev%dstim'%lag]] = 0

    # recombine longer lags
    mydata['prev7_10resp'] = mydata['prev7resp'] + mydata['prev8resp'] + mydata['prev9resp'] + mydata['prev10resp'] 
    mydata['prev7_10stim'] = mydata['prev7stim'] + mydata['prev8stim'] + mydata['prev9stim'] + mydata['prev10stim']
    mydata['prev11_15resp'] = mydata['prev11resp'] + mydata['prev12resp'] + mydata['prev13resp'] + mydata['prev14resp'] + mydata['prev15resp']
    mydata['prev11_15stim'] = mydata['prev11stim'] + mydata['prev12stim'] + mydata['prev13stim'] + mydata['prev14stim'] + mydata['prev15stim']

    return mydata

# ============================================ #
# define the function that will do the work
# ============================================ #

def make_model(mypath, mydata, model_name, trace_id):

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename

    # ============================================ #
    # NO HISTORY FOR MODEL COMPARISON
    # ============================================ #

    if model_name == 'stimcoding_nohist':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'])

    # ============================================ #
    # STIMCODING PREVRESP
    # ============================================ #

    elif model_name == 'stimcoding_dc_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp']})

    elif model_name == 'stimcoding_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'z':['prevresp']})

    # this is the model that will generate most of the figures
    elif model_name == 'stimcoding_dc_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob'],
                'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp'], 'z':['prevresp']})

    # ============================================ #
    # STIMCODING PREVSTIM
    # ============================================ #

    elif model_name == 'stimcoding_dc_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevstim']})

    elif model_name == 'stimcoding_z_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'z':['prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'z':['prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'z':['prevstim']})

    # this is the model that will generate most of the figures
    elif model_name == 'stimcoding_dc_z_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevstim', 'transitionprob'],
                'z':['prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevstim'], 'z':['prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevstim'], 'z':['prevstim']})

    # ============================================ #
    # different non-decision time per coherence
    # ============================================ #

    elif model_name == 'stimcoding_nohist_stcoh':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                                    drift_criterion=True, bias=True, p_outlier=0.05,
                                    include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                                    depends_on={'v': ['coherence'], 't': ['coherence']})

    elif model_name == 'stimcoding_z_prevresp_stcoh':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'z':['prevresp'], 't': ['coherence']})

    elif model_name == 'stimcoding_dc_prevresp_stcoh':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 't': ['coherence']})

    elif model_name == 'stimcoding_dc_z_prevresp_stcoh':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 't': ['coherence'],
                            'dc':['prevresp'], 'z':['prevresp']})

    # ============================================ #
    # also estimate across-trial variability in nondecision time
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevresp_st':

        # include across-trial variability in starting point
        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz', 'st'), group_only_nodes=['sv', 'sz', 'st'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob'],
                'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz', 'st'), group_only_nodes=['sv', 'sz', 'st'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz', 'st'), group_only_nodes=['sv', 'sz', 'st'],
                depends_on={'dc':['prevresp'], 'z':['prevresp']})

    # ============================================ #
    # PHARMA
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevresp_pharma':

        if 'drug' in mydata.columns:
            # get the right variable coding
            mydata = recode_4stimcoding(mydata)

            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'drug'], 'z':['prevresp', 'drug']})
        else:
            m = [] # don't return a model, won't run

    # ============================================ #
    # DOES DRIFT RATE VARIABILITY REDUCE? GROUP
    # ============================================ #

    elif model_name == 'stimcoding_nohist_svgroup':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # add a group to indicate bias magnitude
        mydata['repeat'] 	= (mydata.response == (mydata.prevresp > 0))
        sjrepetition 		= mydata.groupby(['subj_idx'])['repeat'].mean().reset_index()
        sjrepetition['biasgroup'] = pd.qcut(np.abs(sjrepetition['repeat'] - 0.5), 3, labels=False)
        mydata2 = pd.merge(mydata, sjrepetition, on='subj_idx', how='inner')

        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sz'],
                depends_on={'v': ['coherence'], 'sv': ['biasgroup']})
        else:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sz'],
                depends_on={'sv': ['biasgroup']})

    elif model_name == 'stimcoding_dc_z_prevresp_svgroup':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # add a group to indicate bias magnitude
        mydata['repeat'] 	= (mydata.response == (mydata.prevresp > 0))
        sjrepetition 		= mydata.groupby(['subj_idx'])['repeat'].mean().reset_index()
        sjrepetition['biasgroup'] = pd.qcut(np.abs(sjrepetition['repeat'] - 0.5), 3, labels=False)
        mydata2 = pd.merge(mydata, sjrepetition, on='subj_idx', how='inner')

        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 'z':['prevresp'],
                'sv': ['biasgroup']})
        else:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sz'],
                depends_on={'dc':['prevresp'], 'z':['prevresp'],
                'sv': ['biasgroup']})

    # ============================================ #
    # SEPARATE FIT FOR REPEATERS AND ALTERNATORS
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevresp_groupsplit':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # add coding for repeaters and alternators
        mydata['repeat']    = (mydata.response == (mydata.prevresp > 0))
        sjrepetition        = mydata.groupby(['subj_idx'])['repeat'].mean().reset_index()
        sjrepetition['group'] = np.sign(sjrepetition['repeat'] - 0.5)
        mydata2 = pd.merge(mydata, sjrepetition, on='subj_idx', how='inner')

        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'group'], 'z':['prevresp', 'group']})
        else:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'group'], 'z':['prevresp', 'group']})

    # ============================================ #
    # STIMCODING PREVRESP + PREVCORRECT
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevcorrect':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob'],
                'z':['prevresp', 'prevcorrect', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})

    # ============================================ #
    # SUBSAMPLE TO HAVE THE SAME NUMBER OF TRIALS
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevcorrect_subsampled':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # SUBSAMPLE THE DATA FOR PREVIOUS CORRECT, SO THAT FOR
        # EACH PARTICIPANT THERE ARE THE SAME NUMBER OF TRIALS FOR PREVIOUS CORRECT AND ERROR
        trialcounts = mydata.groupby(['subj_idx', 'prevcorrect'])['response'].count().reset_index()
        trialcounts_error = trialcounts.groupby(['subj_idx'])['response'].min().reset_index()

        for name, grouped in mydata.groupby(['subj_idx', 'prevcorrect']):

            # find the lowest number of trials for this subject
            num_trials = trialcounts_error.loc[trialcounts_error['subj_idx'] == name[0], 'response'].item()
            subsampled = grouped.sample(n=num_trials)

            # subsample that many
            if not 'mydata2' in locals():
                mydata2 = subsampled.copy()
            else:
                mydata2 = mydata2.append(subsampled.copy(), ignore_index=True)

        # check that now, the trials are equally numbered
        trialcounts_subsampled = mydata2.groupby(['subj_idx', 'prevcorrect'])['response'].count().reset_index()
        print('trialcounts, before subsampling:')
        print(trialcounts)
        print('trialcounts, after subsampling:')
        print(trialcounts_subsampled)

        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata2, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})


    elif model_name == 'stimcoding_dc_prevcorrect':

        # only let drift bias vary with previous response and previous outcome
        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'prevcorrect']})

    # ============================================ #
    # STIMCODING PREVRESP + PREVCORRECT
    # add a model where previous reward changes drift rate / boundary separation
    # ============================================ #

    elif model_name == 'stimcoding_prevcorrect':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        if len(mydata.coherence.unique()) > 1:
         m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
             drift_criterion=True, bias=True, p_outlier=0.05,
             include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
             depends_on={'v': ['coherence', 'prevcorrect'], 'a':['prevcorrect']})
        else:
         m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
             drift_criterion=True, bias=True, p_outlier=0.05,
             include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
             depends_on={'v': ['prevcorrect'], 'a':['prevcorrect']})

    elif model_name == 'stimcoding_dc_z_PES':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # allow for both choice history effects as well as post-error slowing,
        # so combine 'stimcoding_prevcorrect' and 'stimcoding_dc_z_prevresp'

        if len(mydata.coherence.unique()) > 1:
         m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
             drift_criterion=True, bias=True, p_outlier=0.05,
             include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
             depends_on={'v': ['coherence', 'prevcorrect'], 'a':['prevcorrect'], 
             'dc':['prevresp'], 'z':['prevresp']})
        else:
         m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
             drift_criterion=True, bias=True, p_outlier=0.05,
             include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
             depends_on={'v': ['prevcorrect'], 'a':['prevcorrect'],
             'dc':['prevresp'], 'z':['prevresp']})

    # ============================================ #
    # REPEAT FOR LAG 2
    # ============================================ #

    # this is the model that will generate most of the figures
    elif model_name == 'stimcoding_dc_z_prev2resp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
          m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
              drift_criterion=True, bias=True, p_outlier=0.05,
              include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
              depends_on={'v': ['coherence'], 'dc':['prev2resp', 'transitionprob'],
              'z':['prev2resp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
          m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
              drift_criterion=True, bias=True, p_outlier=0.05,
              include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
              depends_on={'v': ['coherence'], 'dc':['prev2resp'], 'z':['prev2resp']})
        else:
          m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
              drift_criterion=True, bias=True, p_outlier=0.05,
              include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
              depends_on={'dc':['prev2resp'], 'z':['prev2resp']})

    # ============================================ #
    # MULTIPLICATIVE DRIFT BIAS
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevresp_multiplicative':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # this only makes sense with multiple coherence levels
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['coherence', 'prevresp'], 'z':['prevresp']})
        else:
            m = []

    elif model_name == 'stimcoding_dc_prevresp_multiplicative':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # this only makes sense with multiple coherence levels
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['coherence', 'prevresp']})
        else:
            m = []

    # ============================================ #
    # SPLIT BY CONGRUENCE BETWEEN PREVIOUS CHOICE AND CURRENT STIMULUS
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prevresp_congruency':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # add a column coding for the previous congruency
        mydata['stimulus_signed'] = np.sign(mydata.stimulus - 0.2)
        # compute a double (not boolean) to indicate congruence
        mydata['congruent'] = (mydata['prevresp'] == mydata['stimulus_signed']) * 1

        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'congruent'], 'z':['prevresp', 'congruent']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
                depends_on={'dc':['prevresp', 'congruent'], 'z':['prevresp', 'congruent']})

    # ============================================ #
    # MEG DATA
    # ============================================ #

    elif model_name == 'regress_nohist':

        # only stimulus dependence
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + C(transitionprob):prevresp', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + C(transitionprob):prevresp', 'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

 
    # ============================================ #
    # REGRESSION MODELS WITH MULTIPLE LAGS
    # only prevresp, does this improve on AIC?
    # ============================================ #

    elif model_name == 'regress_dc_prevresp_lag1':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_prevresp_lag1':

        z_reg = {'model': 'z ~ 1  + prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, [z_reg, v_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_prevresp_lag1':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        z_reg = {'model': 'z ~ 1  + prevresp ', 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)


    # ============================================ #
    # REGRESSION MODELS WITH MULTIPLE LAGS
    # PREVRESP + PREVSTIM
    # ============================================ #

    elif model_name == 'regress_dc_lag1':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag1':

        z_reg = {'model': 'z ~ 1  + prevresp + prevstim', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, [z_reg, v_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag1':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim', 'link_func': lambda x:x}
        z_reg = {'model': 'z ~ 1  + prevresp + prevstim', 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)
    
    # ============================================ #
    # LAG 2 AND 3
    # ============================================ #

    elif model_name == 'regress_dc_lag2':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag2':

        z_reg = {'model': 'z ~ 1 + prevresp + prevstim + prev2resp + prev2stim', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag2':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim', 'link_func': lambda x:x}
        z_reg = {'model': 'z ~ 1 + prevresp + prevstim + prev2resp + prev2stim', 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 3
    elif model_name == 'regress_dc_lag3':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + prev3resp + prev3stim', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag3':

        z_reg = {'model': 'z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + prev3resp + prev3stim', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag3':

        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + prev3resp + prev3stim', 'link_func': lambda x:x}
        z_reg = {'model': 'z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + prev3resp + prev3stim', 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 4
    elif model_name == 'regress_dc_lag4':

        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag4':

        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag4':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 5
    elif model_name == 'regress_dc_lag5':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag5':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag5':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 6
    elif model_name == 'regress_dc_lag6':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag6':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag6':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 7
    elif model_name == 'regress_dc_lag7':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag7':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag7':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # RECODE THE MODEL!!!!
    elif model_name == 'regress_dcz_lag7_recode':

        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prev1resp_correct + prev1resp_error 
        + prev2resp_correct + prev2resp_error + prev3resp_correct + prev3resp_error
        + prev4resp_correct + prev4resp_error + prev5resp_correct + prev5resp_error
        + prev6resp_correct + prev6resp_error + prev7resp_correct + prev7resp_error""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prev1resp_correct + prev1resp_error 
        + prev2resp_correct + prev2resp_error + prev3resp_correct + prev3resp_error
        + prev4resp_correct + prev4resp_error + prev5resp_correct + prev5resp_error
        + prev6resp_correct + prev6resp_error + prev7resp_correct + prev7resp_error""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # LAG 8-11
    elif model_name == 'regress_dc_lag7-10':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7_10resp + prev7_10stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag7-10':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7_10resp + prev7_10stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag7-10':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7_10resp + prev7_10stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7_10resp + prev7_10stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)


    # LAG 11_15
    elif model_name == 'regress_dc_lag11-15':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7_10resp + prev7_10stim 
        + prev11_15resp + prev11_15stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag11-15':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7_10resp + prev7_10stim
         + prev11_15resp + prev11_15stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag11-15':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7resp + prev7stim + prev7_10resp + prev7_10stim
         + prev11_15resp + prev11_15stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7_10resp + prev7_10stim
        + prev11_15resp + prev11_15stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # ============================================ #
    # INCLUDE TRIAL-BY-TRIAL MEG REGRESSORS
    # ============================================ #


   # elif model_name == 'regress_dc_z_visualgamma':

   #      z_reg = {'model': 'z ~ 1 + visualgamma', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + visualgamma', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

   #  elif model_name == 'regress_dc_z_motorslope':

   #      z_reg = {'model': 'z ~ 1 + motorslope', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + motorslope', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

   #  elif model_name == 'regress_dc_z_motorstart':

   #      z_reg = {'model': 'z ~ 1 + motorbeta', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + motorbeta', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)


   #  elif model_name == 'regress_dc_z_prevresp_visualgamma':

   #      z_reg = {'model': 'z ~ 1 + prevresp + visualgamma', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + prevresp + visualgamma', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

   #  elif model_name == 'regress_dc_z_prevresp_motorslope':

   #      z_reg = {'model': 'z ~ 1 + prevresp + motorslope', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + prevresp + motorslope', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

   #  elif model_name == 'regress_dc_z_prevresp_motorstart':

   #      z_reg = {'model': 'z ~ 1 + prevresp + motorbeta', 'link_func': z_link_func}
   #      v_reg = {'model': 'v ~ 1 + stimulus + prevresp + motorbeta', 'link_func': lambda x:x}
   #      reg_both = [z_reg, v_reg]

   #      m = hddm.HDDMRegressor(mydata, reg_both,
   #      include=['z', 'sv'], group_only_nodes=['sv'],
   #      group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

    # ============================================ #
    # END OF FUNCTION THAT CREATES THE MODEL
    # ============================================ #

    return m
