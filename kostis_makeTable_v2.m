function kostis_makeTable_v2

close all; clc;
kostisPath = '~/Data/HDDM/Anke_MEG_transition/KostisFits';

% GRAB MOTION ENERGY TO DETERMINE NR OF TRIALS USED FOR THE FIT
load(sprintf('%s/motionEnergyData_AnkeMEG.mat', kostisPath));
ss=unique(data.behavior.subj_idx);
for s=1:length(ss);
    indx=find(data.behavior.subj_idx==ss(s) & data.behavior.RT>0.25 & data.behavior.coherence~=81  &  ~isnan(data.behavior.prevresp));
    notrials(s)=length(indx);
end

% function for BIC computation - treating all subjects as a super subject
% 2*sum(LL(:,1))+(length(ss)*5).*log(sum(notrials))
ll2bic = @(ll, p, n) repmat((2*sum(ll) + length(n)*p.*log(sum(n))), length(n), 1);

% ========================================== %
% 1. DDM values
% ========================================== %

% params: naive model
% params2: DC
% params3: starting
% params4: ramping 1-par
% params5: ramping 2-pars
% params6: DC+starting

% in all models the parameters are:
% column 1: boundary: multiply by 10
% column 2: scale controlling signal-to-noise (noise fixed at std=1):
% multiply by 60

% column 3: T0 expressed in seconds: multiply by 0.5 seconds
% column 4: drift-rate variability: multiply by 3
% column 5: starting point variability: multiply by 0.75. Parameter
% expresses % relative to threshold i.e. 0.75 means 75% of the value of the
% trheshold
% for params6: column 6: drift-criterion, column 7: starting point

% likelihoods for each subject
% rows are: 1-naive, 2-DC, 3-startint point
% 4-ramping with slope, 5-ramping with 2 params
% 6- both DC and starting point;

clearvars -except kostisPath notrials ll2bic 
load(sprintf('%s/allmodels.mat', kostisPath));

params_ddm          = array2table([params ll2bic(outgf(:, 1), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'bic'});
params_ddm_dc       = array2table([params2 ll2bic(outgf(:, 2), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'dcbias', 'bic'});
params_ddm_sp       = array2table([params3 ll2bic(outgf(:, 3), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'zbias', 'bic'});
params_ddm_sp_dc    = array2table([params6 ll2bic(outgf(:, 6), 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'dcbias', 'zbias', 'bic'});
%column 6: in params2, params3, params4 is the biasing parameter.
% multiply by 1,1 and 5 respectively for those models.
params_ddm_rp       = array2table([params4 ll2bic(outgf(:, 4), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'slope', 'bic'});
params_ddm_rp.slope = params_ddm_rp.slope * 5;

%for params5, column 6: slope (multiply by sign of offset), column 7: offset
params_ddm_rp2      = array2table([params5 ll2bic(outgf(:, 5), 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'slope', 'offset', 'bic'});
params_ddm_rp2.slope = params_ddm_rp2.slope .* sign(params_ddm_rp2.offset);

params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ddmK_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_z_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_dc.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_dc_', x)), params_ddm_dc.Properties.VariableNames, 'un', 0);
params_ddm_sp_dc.Properties.VariableNames   = cellfun((@(x) cat(2, 'ddmK_dcz_', x)), params_ddm_sp_dc.Properties.VariableNames, 'un', 0);
params_ddm_rp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_rp_', x)), params_ddm_rp.Properties.VariableNames, 'un', 0);
params_ddm_rp2.Properties.VariableNames     = cellfun((@(x) cat(2, 'ddmK_rp2_', x)), params_ddm_rp2.Properties.VariableNames, 'un', 0);

alltables{1} = cat(2, params_ddm, params_ddm_sp, params_ddm_dc, params_ddm_sp_dc, params_ddm_rp, params_ddm_rp2);

% ========================================== %
% 2. DDM WITH collapsing bounds
% ========================================== %

% load DDMCol_allmodels.mat % collapsing DDM
% [a1,x1]=corr(kk(:,2),params2(:,end),'Type','Spearman');% drift-criterion
% [a2,x2]=corr(kk(:,2),params3(:,end),'Type','Spearman'); % starting point
% [a3a,x3a]=corr(kk(:,2),params4(:,end),'Type','Spearman');% starting point in hybrid
% [a3b,x3b]=corr(kk(:,2),params4(:,end-1),'Type','Spearman'); % drift criterion hybrid

clearvars -except kostisPath notrials ll2bic alltables
load(sprintf('%s/DDMCol_allmodels.mat', kostisPath));

% assume that the size of individual bound collapse (hyperbolic) is in the
% same column in the parameter matrices as the drift rate variability in
% the vanilla DDMs above?
params_ddm          = array2table([params ll2bic(outgf(:, 1), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'boundcollapse', 'bsp', 'bic'});
params_ddm_dc       = array2table([params2 ll2bic(outgf(:, 2), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'boundcollapse', 'bsp', 'dcbias', 'bic'});
params_ddm_sp       = array2table([params3 ll2bic(outgf(:, 3), 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'boundcollapse', 'bsp', 'zbias', 'bic'});
params_ddm_sp_dc    = array2table([params4 ll2bic(outgf(:, 4), 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'boundcollapse', 'bsp', 'dcbias', 'zbias', 'bic'});

params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ddmColl_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmColl_z_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_dc.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmColl_dc_', x)), params_ddm_dc.Properties.VariableNames, 'un', 0);
params_ddm_sp_dc.Properties.VariableNames   = cellfun((@(x) cat(2, 'ddmColl_dcz_', x)), params_ddm_sp_dc.Properties.VariableNames, 'un', 0);

alltables{2} = cat(2, params_ddm, params_ddm_sp, params_ddm_dc, params_ddm_sp_dc);

% ========================================== %
% 3. then take the O-U values
% ========================================== %

% params: naive model
% params2: Inout bias
% params3: OU-asymmetry
% params4: starting point

% in all models the parameters are:
% column 1: boundary: multiply by 30
% column 2: scale controlling signal-to-noise (noise fixed at std=1):
% multiply by 100

% column 3: T0 expressed in seconds: multiply by 0.7 seconds
% column 4: OU parameter: multiply by 15
% column 5: starting point variability: multiply by 0.75. Parameter
% expresses % relative to threshold i.e. 0.75 means 75% of the value of the
% trheshold

clearvars -except kostisPath notrials ll2bic alltables
load(sprintf('%s/OUallmodels.mat', kostisPath));

params_ddm          = array2table([params ll2bic(outgf(:, 1), 5, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'bic'});
params_ddm_input    = array2table([params2 ll2bic(outgf(:, 2), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'inputbias', 'bic'});
params_ddm_lambda   = array2table([params3 ll2bic(outgf(:, 3), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'lambdabias', 'bic'});
params_ddm_sp       = array2table([params4 ll2bic(outgf(:, 4), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'spbias', 'bic'});

% column 6: in params2, params3, params4 is the biasing parameter.
% multiply by 5,1 and 5 respectively for those models.
params_ddm_input.inputbias = params_ddm_input.inputbias * 5;
params_ddm_sp.spbias = params_ddm_sp.spbias * 5;

% put into one big table
params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ouK_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ouK_sp_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_input.Properties.VariableNames   = cellfun((@(x) cat(2, 'ouK_input_', x)), params_ddm_input.Properties.VariableNames, 'un', 0);
params_ddm_lambda.Properties.VariableNames  = cellfun((@(x) cat(2, 'ouK_lambda_', x)), params_ddm_lambda.Properties.VariableNames, 'un', 0);

alltables{3} = cat(2, params_ddm, params_ddm_sp, params_ddm_input, params_ddm_lambda);

% ========================================== %
% 4. O-U, accumulation during stimulus
% ========================================== %

% params: naive model
% params2: Inout bias
% params3: OU-asymmetry
% params4: starting point

% in all models the parameters are:
% column 1: boundary: multiply by 30
% column 2: scale controlling signal-to-noise (noise fixed at std=1):
% multiply by 100

% column 3: T0 expressed in seconds: multiply by 0.7 seconds
% column 4: OU parameter: multiply by 15
% column 5: starting point variability: multiply by 0.75. Parameter
% expresses % relative to threshold i.e. 0.75 means 75% of the value of the
% trheshold

clearvars -except kostisPath notrials ll2bic alltables
load(sprintf('%s/OUDallmodels.mat', kostisPath));

params_ddm          = array2table([params ll2bic(outgf(:, 1), 5, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'bic'});
params_ddm_input    = array2table([params2 ll2bic(outgf(:, 2), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'inputbias', 'bic'});
params_ddm_lambda   = array2table([params3 ll2bic(outgf(:, 3), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'lambdabias', 'bic'});
params_ddm_sp       = array2table([params4 ll2bic(outgf(:, 4), 6, notrials)], 'variablenames', {'boundary', 'scale', 'T0', 'lambda', 'bsp', 'spbias', 'bic'});

% column 6: in params2, params3, params4 is the biasing parameter.
% multiply by 5,1 and 5 respectively for those models.
params_ddm_input.inputbias = params_ddm_input.inputbias * 5;
params_ddm_sp.spbias = params_ddm_sp.spbias * 5;

% put into one big table
params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ouD_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ouD_sp_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_input.Properties.VariableNames   = cellfun((@(x) cat(2, 'ouD_input_', x)), params_ddm_input.Properties.VariableNames, 'un', 0);
params_ddm_lambda.Properties.VariableNames  = cellfun((@(x) cat(2, 'ouD_lambda_', x)), params_ddm_lambda.Properties.VariableNames, 'un', 0);

alltables{4} = cat(2, params_ddm, params_ddm_sp, params_ddm_input, params_ddm_lambda);

% ========================================== %
% CONCATENATE ALL INTO A TABLE
% ========================================== %

% add repetition, computed from all of Anke's trials
load(sprintf('%s/history.mat', kostisPath));
kk(6, :) = []; % remove this subject
alltables{end}.repetitionK = kk(:, 2); % add to the big table

% load the main results file
results = readtable('~/Data/HDDM/summary/Anke_MEG_transition/allindividualresults.csv');
results = results(results.session == 0, :);
results = cat(2, results, alltables{:});
writetable(results, '~/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');

corrplot(results, {'repetitionK', 'ddmK_z_zbias', 'ddmK_dc_dcbias', 'ddmK_dcz_zbias', 'ddmK_dcz_dcbias'});
print(gcf, '-dpdf', '~/Data/serialHDDM/kostisData_overview_DDM.pdf');

close all;
corrplot(results, {'repetitionK', 'ddmK_rp_slope', 'ddmK_dc_dcbias', ...
    'ddmK_rp2_slope', 'ddmK_rp2_offset'});
print(gcf, '-dpdf', '~/Data/serialHDDM/kostisData_overview_DDM_ramp.pdf');

close all;
corrplot(results, {'repetitionK', 'repetition'})
print(gcf, '-dpdf', '~/Data/serialHDDM/repetition_comparison.pdf');


end