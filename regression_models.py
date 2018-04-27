    # ============================================ #
    # REGRESSION MODULATION
    # Nienborg @ SfN: ((s*v) + vbias)dt / (s*(v+vbias))dt, does it matter and what are the differential predictions?
    # ============================================ #

    elif model_name == 'regress_dc_prevresp':

        if 'transitionprob' in mydata.columns:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp:C(transitionprob) ', 'link_func': lambda x:x}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    elif model_name == 'regress_dc2_prevresp':

        v_reg = {'model': 'v ~ (1 + prevresp)*stimulus', 'link_func': lambda x:x}
        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # ============================================ #
    # MODULATION BY PUPIL AND RT
    # ============================================ #

    elif model_name == 'regress_dc_prevresp_prevrt':

        if 'transitionprob' in mydata.columns:
         v_reg = {'model': 'v ~ 1 + stimulus + C(transitionprob):prevresp + ' \
         'C(transitionprob):prevrt + C(transitionprob):prevresp:prevrt', 'link_func': lambda x:x}
        else:
         v_reg = {'model': 'v ~ 1 + stimulus + prevresp*prevrt', 'link_func': lambda x:x}

        m = hddm.HDDMRegressor(mydata, v_reg,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

    # ================================================== #
    # THE MOST IMPORTANT MODEL FOR THE UNCERTAINTY MODULATION
    # ================================================== #

    elif model_name == 'regress_dc_z_prevresp_prevrt':

        mydata = mydata.dropna(subset=['prevresp','prevrt'])

        # subselect data
        if 'transitionprob' in mydata.columns:
             v_reg = {'model': 'v ~ 1 + stimulus + C(transitionprob):prevresp + ' \
             'C(transitionprob):prevrt + C(transitionprob):prevresp:prevrt', 'link_func': lambda x:x}
             z_reg = {'model': 'z ~ 1 + prevresp:C(transitionprob) +' \
             'prevresp:prevrt:C(transitionprob) + prevrt:C(transitionprob)',
             'link_func': z_link_func}
        else:
            v_reg = {'model': 'v ~ 1 + stimulus + prevresp*prevrt',
            'link_func': lambda x:x}
            z_reg = {'model': 'z ~ 1 + prevresp*prevrt',
            'link_func': z_link_func}

        reg_both = [v_reg, z_reg]
        m = hddm.HDDMRegressor(mydata, reg_both,
        include=['z', 'sv'], group_only_nodes=['sv'],
        group_only_regressors=False, keep_regressor_trace=False, p_outlier=0.05)

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
        group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)
        
        
        
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
            group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

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
            group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

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
            group_only_regressors=False, keep_regressor_trace=False,  p_outlier=0.05)

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
        
