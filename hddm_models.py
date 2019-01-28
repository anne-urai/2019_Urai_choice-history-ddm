#!/usr/bin/env python
# encoding: utf-8

"""
Code to fit the history-dependent drift diffusion models described in
Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595

MIT License
Copyright (c) Anne Urai, 2018
anne.urai@gmail.com
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
    mydata['prev7-10resp'] = mydata['prev7resp'] + mydata['prev8resp'] + mydata['prev9resp'] + mydata['prev10resp'] 
    mydata['prev7-10stim'] = mydata['prev7stim'] + mydata['prev8stim'] + mydata['prev9stim'] + mydata['prev10stim']
    mydata['prev11-15resp'] = mydata['prev11resp'] + mydata['prev12resp'] + mydata['prev13resp'] + mydata['prev14resp'] + mydata['prev15resp']
    mydata['prev11-15stim'] = mydata['prev11stim'] + mydata['prev12stim'] + mydata['prev13stim'] + mydata['prev14stim'] + mydata['prev15stim']

    shell()
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

    # this is the model that will generate most of the figures
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

    elif model_name == 'regress_dc_z_visualgamma':

        z_reg = {'model': 'z ~ 1 + visualgamma', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + visualgamma', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_motorslope':

        z_reg = {'model': 'z ~ 1 + motorslope', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + motorslope', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_motorstart':

        z_reg = {'model': 'z ~ 1 + motorbeta', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + motorbeta', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)


    elif model_name == 'regress_dc_z_prevresp_visualgamma':

        z_reg = {'model': 'z ~ 1 + prevresp + visualgamma', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + visualgamma', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp_motorslope':

        z_reg = {'model': 'z ~ 1 + prevresp + motorslope', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + motorslope', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp_motorstart':

        z_reg = {'model': 'z ~ 1 + prevresp + motorbeta', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus + prevresp + motorbeta', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)


    # ============================================ #
    # REGRESSION MODELS WITH MULTIPLE LAGS
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

        mydata['prevresp_correct']

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7resp + prev7stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)


    # LAG 8-11
    elif model_name == 'regress_dc_lag7-10':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7-10resp + prev7-10stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag7-10':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7-10resp + prev7-10stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag7-10':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7-10resp + prev7-10stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7-10resp + prev7-10stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)


    # LAG 12-15
    elif model_name == 'regress_dc_lag11-15':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7-10resp + prev7-10stim 
        + prev11-15resp + prev11-15stim""", 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_z_lag11-15':
        mydata = add_morelags(mydata)

        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7-10resp + prev7-10stim
         + prev11-15resp + prev11-15stim""", 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dcz_lag12-15':
        mydata = add_morelags(mydata)

        v_reg = {'model': """v ~ 1 + stimulus + prevresp + prevstim + prev2resp + prev2stim +
         prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
         + prev6resp + prev6stim + prev7resp + prev7stim + prev7-10resp + prev7-10stim
         + prev11-15resp + prev11-15stim""", 'link_func': lambda x:x}
        z_reg = {'model': """z ~ 1 + prevresp + prevstim + prev2resp + prev2stim + 
        prev3resp + prev3stim + prev4resp + prev4stim + prev5resp + prev5stim
        + prev6resp + prev6stim + prev7-10resp + prev7-10stim
        + prev11-15resp + prev11-15stim""", 'link_func': z_link_func}
        m = hddm.HDDMRegressor(mydata, [v_reg, z_reg],
                               include=['z', 'sv'], group_only_nodes=['sv'],
                               group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # ============================================ #
    # END OF FUNCTION THAT CREATES THE MODEL
    # ============================================ #

    return m
