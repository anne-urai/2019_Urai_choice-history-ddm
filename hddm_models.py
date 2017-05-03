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

# ============================================ #
# define the function that will do the work
# ============================================ #

def make_model(mypath, model_name, trace_id):

    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    print model_filename

    # get the csv file for this dataset
    filename    = fnmatch.filter(os.listdir(mypath), '*.csv')
    mydata      = hddm.load_csv(os.path.join(mypath, filename[0]))

    # prepare link function for the regression models
    def z_link_func(x, data=mydata):
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
        mydata.stimulus  = np.sign(mydata.stimulus)
        # for stimcoding, the two identities should be 0 and 1
        mydata.ix[mydata['stimulus']==-1,'stimulus'] = 0
        return mydata

    # ============================================ #
    # NO HISTORY FOR MODEL COMPARISON
    # ============================================ #

    if model_name == 'stimcoding_nohist':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'])

    if model_name == 'regress_nohist':

        # only stimulus dependence
        v_reg = {'model': 'v ~ 1 + stimulus', 'link_func': lambda x:x}

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # STIMCODING PREVRESP
    # ============================================ #

    if model_name == 'stimcoding_dc_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp']})

    if model_name == 'stimcoding_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'z':['prevresp']})

    if model_name == 'stimcoding_dc_z_prevresp':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'transitionprob'],
                'z':['prevresp', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp'], 'z':['prevresp']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp'], 'z':['prevresp']})

    # ============================================ #
    # STIMCODING PREVRESP + PREVSTIM
    # ============================================ #

    if model_name == 'stimcoding_dc_prevresp_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability and include coherence-dependence of drift rate
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            print "splitting by coherence"
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevstim']})

    if model_name == 'stimcoding_z_prevresp_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'z':['prevresp', 'prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'z':['prevresp', 'prevstim']})

    if model_name == 'stimcoding_dc_z_prevresp_prevstim':

        # get the right variable coding
        mydata = recode_4stimcoding(mydata)

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevstim', 'transitionprob'],
                'z':['prevresp', 'prevstim', 'transitionprob']})
        elif len(mydata.coherence.unique()) > 1: # Ankes neutral condition
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'v': ['coherence'], 'dc':['prevresp', 'prevstim'], 'z':['prevresp', 'prevstim']})
        else:
            m = hddm.HDDMStimCoding(mydata, stim_col='stimulus', split_param='v',
                drift_criterion=True, bias=True, p_outlier=0.05,
                include=('sv'), group_only_nodes=['sv'],
                depends_on={'dc':['prevresp', 'prevstim'], 'z':['prevresp', 'prevstim']})

    # ============================================ #
    # REGRESS PREVRESP
    # ============================================ #

    if model_name == 'regress_dc_prevresp':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_z_prevresp':

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
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_z_prevresp':

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
        group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # REGRESS PREVRESP + PREVSTIM
    # ============================================ #

    if model_name == 'regress_dc_prevresp_prevstim':

        # for Anke's data, also split by transition probability
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob) +' \
                'prevstim:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_z_prevresp_prevstim':

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
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_z_prevresp_prevstim':

        if 'transitionprob' in mydata.columns:
            z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) + prevstim:C(transitionprob)',
                'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob) + prevstim:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            z_reg = {'model': 'z ~ 1 + prevresp + prevstim', 'link_func': z_link_func}
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp + prevstim', 'link_func': lambda x:x}
        reg_both = [z_reg, v_reg]

        # specify that we want individual parameters for all regressors, see email Gilles 22.02.2017
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # SESSION DEPENDENCE
    # ============================================ #

    if model_name == 'regress_dc_prevresp_prevstim_vasessions':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) +' \
                'prevstim:C(transitionprob)', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim', 'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrespsessions':

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
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevstim_vasessions_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + ' \
                'prevresp:C(transitionprob) + prevstim:C(transitionprob) + ' \
                'prevresp:prevpupil:C(transitionprob) + prevstim:prevpupil:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim + ' \
                'prevresp:prevpupil + prevstim:prevpupil',
                'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrt':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) +' \
                'prevstim:C(transitionprob) + prevresp:prevrt:C(transitionprob) + prevstim:prevrt:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prevstim + prevresp:prevrt + prevstim:prevrt',
                'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prevresp_prevstim_vasessions_prevrt_prevpupil':

        # subselect data
        mydata = mydata.dropna(subset=['prevpupil'])
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
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
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x}
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    # ============================================ #
    # INCLUDE LAGS FURTHER INTO THE PAST
    # ============================================ #

    if model_name == 'regress_dc_prev2resp_prev2stim':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) +' \
                'prevstim:C(transitionprob) + prev2resp:C(transitionprob) + prev2stim:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prev2resp + prevstim + prev2stim', 'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    if model_name == 'regress_dc_prev3resp_prev3stim':

        # subselect data
        mydata = balance_designmatrix(mydata)

        # boundary separation and drift rate will change over sessions
        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp:C(transitionprob) +' \
                'prevstim:C(transitionprob) + prev2resp:C(transitionprob) + prev2stim:C(transitionprob) + ' \
                'prev3resp:C(transitionprob) + prev3stim:C(transitionprob)',
                'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus:C(session) + prevresp + prev2resp + ' \
            'prevstim + prev2stim + prev3resp + prev3stim', 'link_func': lambda x:x}
        a_reg = {'model': 'a ~ 1 + C(session)', 'link_func': lambda x:x} # boundary separation as a function of sessions
        reg_both = [v_reg, a_reg]

        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, p_outlier=0.05)

    # END OF FUNCTION THAT CREATES THE MODEL
    return m
