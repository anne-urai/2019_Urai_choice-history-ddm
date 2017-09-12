#!/usr/bin/env python
# encoding: utf-8

"""
Anne Urai, 2017
takes input arguments from stopos
Important: on Cartesius, call module load python/2.7.9 before running
(the only environment where HDDM is installed)
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

    # also create a binary prevcorrect
    mydata['prev2correct']     = mydata.prevresp
    mydata.prev2correct[mydata.prev2resp != mydata.prev2stim] = 0
    mydata.prev2correct[mydata.prev2resp == mydata.prev2stim] = 1
    # also create a binary prevcorrect
    mydata['prev3correct']     = mydata.prevresp
    mydata.prev3correct[mydata.prev3resp != mydata.prev3stim] = 0
    mydata.prev3correct[mydata.prev3resp == mydata.prev3stim] = 1

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
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'])

    if model_name == 'stimcoding_nohist_onlyz':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'])

    elif model_name == 'stimcoding_nohist_onlydc':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'])

    elif model_name == 'regress_nohist':

        # only stimulus dependence
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg,
            include=['z', 'sv'], group_only_nodes=['sv'],
            group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # STIMCODING PREVRESP
    # ============================================ #

    elif model_name == 'stimcoding_dc_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp']})

    elif model_name == 'stimcoding_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'z':['prevresp']})

    elif model_name == 'stimcoding_dc_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob'],
                'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp'], 'z':['prevresp']})

    elif model_name == 'stimcoding_dc_z_prevresp_pharma':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
            drift_criterion=True, bias=True, p_outlier=0.05,
            include=('sv'), group_only_nodes=['sv'],
            depends_on={'dc':['prevresp', 'drug'], 'z':['prevresp', 'drug']})

    # ============================================ #
    # STIMCODING PREVRESP + PREVCORRECT
    # ============================================ #

    elif model_name == 'stimcoding_dc_prevcorrect':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            print "splitting by coherence"
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevcorrect']})

    elif model_name == 'stimcoding_z_prevcorrect':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'prevcorrect', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'z':['prevresp', 'prevcorrect']})

    elif model_name == 'stimcoding_dc_z_prevcorrect':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob'],
                'z':['prevresp', 'prevcorrect', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevcorrect'], 'z':['prevresp', 'prevcorrect']})

    # ============================================ #
    # REGRESSION MODULATION
    # ============================================ #

    elif model_name == 'regress_dc_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)',
            'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # MODULATION BY PUPIL AND RT
    # ============================================ #

    elif model_name == 'regress_dc_z_prevresp_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus+ ' \
                'prevresp:C(transitionprob)  + ' \
                'prevresp:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) +' \
                'prevresp:prevrt:C(transitionprob) ',
                'link_func': z_link_func}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevresp:prevrt',
            'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp  +' \
            'prevresp:prevrt',
            'link_func': z_link_func}
        reg_both = [v_reg, z_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp_prevstim_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
            'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
            'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob) + ' \
            'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
            'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prevstim:C(transitionprob) +' \
            'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob) +' \
            'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
            'link_func': z_link_func}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + ' \
            'prevresp:prevrt + prevstim:prevrt + prevresp:prevpupil + prevstim:prevpupil',
            'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp + prevstim +' \
            'prevresp:prevpupil + prevstim:prevpupil +' \
            'prevresp:prevrt + prevstim:prevrt',
            'link_func': z_link_func}
        reg_both = [v_reg, z_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # MULTIPLE LAGS
    # ============================================ #

    elif model_name == 'regress_dc_z_prev2resp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)+ prev2resp:C(transitionprob)',
            'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)  + prev2resp:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp  + prev2resp ', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prev2resp ', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # subselect data
        mydata = mydata.dropna(subset=['prev2resp'])

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_prev3resp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prev2resp:C(transitionprob) + prev3resp:C(transitionprob)',
            'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob) + prev2resp:C(transitionprob) + prev3resp:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp + prev2resp + prev3resp ', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prev2resp + prev3resp ', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        mydata = mydata.dropna(subset=['prev2resp', 'prev3resp'])

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # DONT USE after this
    # ============================================ #

    elif model_name == 'stimcoding_dc_z_prev2correct':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob', 'prev2resp', 'prev2correct'],
                'z':['prevresp', 'prevcorrect', 'transitionprob','prev2resp', 'prev2correct']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct'], 'z':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct'], 'z':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct']})

    elif model_name == 'stimcoding_dc_z_prev3correct':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'transitionprob', 'prev2resp', 'prev2correct', 'prev3resp', 'prev3correct'],
                'z':['prevresp', 'prevcorrect', 'transitionprob','prev2resp', 'prev2correct', 'prev3resp', 'prev3correct']})
        elif len(mydata.coherence.unique()) > 1:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct', 'prev3resp', 'prev3correct'], 'z':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct', 'prev3resp', 'prev3correct']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct'], 'z':['prevresp', 'prevcorrect', 'prev2resp', 'prev2correct']})

    # ============================================ #
    # REGRESS PREVRESP
    # ============================================ #

    elif model_name == 'regress_dc_prevresp':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True, p_outlier=0.05)

    elif model_name == 'regress_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)',
            'link_func': z_link_func}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True, p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob)',
            'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # REGRESS PREVRESP + PREVSTIM
    # ============================================ #

    elif model_name == 'regress_dc_prevresp_prevstim':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob) +' \
            'prevstim:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_z_prevresp_prevstim':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prevstim:C(transitionprob)',
            'link_func': z_link_func}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp + prevstim', 'link_func': z_link_func}
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # RT AND PUPIL MODULATION onto DC
    # ============================================ #

    elif model_name == 'regress_dc_prevresp_prevstim_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})
            mydata = balance_designmatrix(mydata)

        # in Anke's data, vary both dc and z
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
                'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + ' \
                'prevresp:prevpupil + prevstim:prevpupil', 'link_func': lambda x:x}
        reg_both = [v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
                'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + prevresp:prevrt + prevstim:prevrt',
            'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})

        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
                'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob) + ' \
                'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + ' \
                'prevresp:prevrt + prevstim:prevrt + prevresp:prevpupil + prevstim:prevpupil',
                'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # LEARNING EFFECTS, V AND A CHANGE WITH SESSION
    # ============================================ #

    elif model_name == 'regress_dc_z_prevresp_prevstim_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})
            mydata = balance_designmatrix(mydata)

        # in Anke's data, vary both dc and z
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
                'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
                'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prevstim:C(transitionprob) +' \
                'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
                'link_func': z_link_func}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim + ' \
                'prevresp:prevpupil + prevstim:prevpupil',
                'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp + prevstim +' \
                'prevresp:prevpupil + prevstim:prevpupil',
                'link_func': z_link_func}
        reg_both = [v_reg, z_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # SESSION DEPENDENCE
    # ============================================ #

    elif model_name == 'regress_dc_prevresp_prevstim_vasessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) + prevstim:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_z_prevresp_prevstim_vasessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) + prevstim:C(transitionprob)', 'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prevstim:C(transitionprob)', 'link_func': z_link_func}
            reg_both = [v_reg, a_reg, z_reg]
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim', 'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp+ prevstim', 'link_func': z_link_func}
        reg_both = [v_reg, a_reg, z_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrespsessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            raise ValueError('Do not fit session-specific serial bias on Anke''s data')
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(session) + prevstim:C(session)', 'link_func': lambda x:x}
            a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_vasessions_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + ' \
            'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
            'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim + ' \
            'prevresp:prevpupil + prevstim:prevpupil',
            'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + ' \
            'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
            'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim + prevresp:prevrt + prevstim:prevrt',
            'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})

        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + ' \
            'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
            'prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob) +' \
            'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
            'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim + ' \
            'prevresp:prevrt + prevstim:prevrt + prevresp:prevpupil + prevstim:prevpupil',
            'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # DIFFERENTIATE PREVIOUS ERROR + CORRECT
    # ============================================ #

    elif model_name == 'regress_dc_prevcorrect':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + preverror:C(transitionprob) +' \
            'prevcorrect:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + preverror + prevcorrect', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevcorrect_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})
            mydata = balance_designmatrix(mydata)

        # in Anke's data, vary both dc and z
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevcorrect:C(transitionprob) + preverror:C(transitionprob) + ' \
                'prevcorrect:prevpupil:C(transitionprob) + preverror:prevpupil:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevcorrect + preverror + ' \
                'prevcorrect:prevpupil + preverror:prevpupil', 'link_func': lambda x:x}
        reg_both = [v_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevcorrect_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevcorrect:C(transitionprob) + preverror:C(transitionprob) + ' \
                'prevcorrect:prevrt:C(transitionprob) + preverror:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevcorrect + preverror + prevcorrect:prevrt + preverror:prevrt',
            'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    elif model_name == 'regress_dc_prevcorrect_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        if len(mydata.session.unique()) < max(mydata.session.unique()):
            mydata["session"] = mydata["session"].map({1:1, 5:2})

        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + ' \
                'prevcorrect:C(transitionprob) + preverror:C(transitionprob) + ' \
                'prevcorrect:prevpupil:C(transitionprob) + preverror:prevpupil:C(transitionprob) + ' \
                'prevcorrect:prevrt:C(transitionprob) + preverror:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevcorrect + preverror + ' \
                'prevcorrect:prevrt + preverror:prevrt + prevcorrect:prevpupil + preverror:prevpupil',
                'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=True,  p_outlier=0.05)

    # ============================================ #
    # END OF FUNCTION THAT CREATES THE MODEL
    # ============================================ #

    return m
