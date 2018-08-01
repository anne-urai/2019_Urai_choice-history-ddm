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

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
            depends_on={'dc':['prevresp', 'drug'], 'z':['prevresp', 'drug']})

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
    # END OF FUNCTION THAT CREATES THE MODEL
    # ============================================ #

    return m

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
    # END OF FUNCTION THAT CREATES THE MODEL
    # ============================================ #

    return m
