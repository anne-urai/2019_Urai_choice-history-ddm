function kostis_makeTable

close all; clc;
global mypath datasets datasetnames
kostisPath = '~/Data/HDDM/Anke_MEG_neutral/materials_K/models';

% GRAB MOTION ENERGY TO DETERMINE NR OF TRIALS USED FOR THE FIT
load(sprintf('%s/motionEnergyData_AnkeMEG.mat', kostisPath));
ss=unique(data.behavior.subj_idx);
for s=1:length(ss);
    indx=find(data.behavior.subj_idx==ss(s) & data.behavior.RT>0.25 & data.behavior.coherence~=81  &  ~isnan(data.behavior.prevresp));
    notrials(s)=length(indx);
end

% function for BIC computation
ll2bic = @(ll, p, n) 2*ll'+p.*log(n');

%% 1. first take the O-U values
load(sprintf('%s/cartesius_OU/DDM.mat', kostisPath));
params_ddm = array2table([params ll2bic(Gf, 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'bic'});

load(sprintf('%s/cartesius_OU/DDM_SP.mat', kostisPath));
params_ddm_sp = array2table([params ll2bic(Gf, 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'spbias', 'bic'});

load(sprintf('%s/cartesius_OU/DDM_DC.mat', kostisPath));
params_ddm_input = array2table([params ll2bic(Gf, 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'inputbias', 'bic'});

load(sprintf('%s/cartesius_OU/DDM_AP.mat', kostisPath));
params_ddm_lambda = array2table([params ll2bic(Gf, 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'lambdabias', 'bic'});

% put into one big table
params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ouK_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ouK_sp_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_input.Properties.VariableNames   = cellfun((@(x) cat(2, 'ouK_input_', x)), params_ddm_input.Properties.VariableNames, 'un', 0);
params_ddm_lambda.Properties.VariableNames  = cellfun((@(x) cat(2, 'ouK_lambda_', x)), params_ddm_lambda.Properties.VariableNames, 'un', 0);

ou_table = cat(2, params_ddm, params_ddm_sp, params_ddm_input, params_ddm_lambda);

% add repetition, computed from all of Anke's trials
load(sprintf('%s/history.mat', kostisPath));
kk(6, :) = []; % remove this subject
ou_table.repetition_alldata = kk(:, 2);

%% 2. DDM values
load(sprintf('%s/DDM.mat', kostisPath));
params_ddm = array2table([params ll2bic(Gf, 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'bic'});

load(sprintf('%s/DDM_SP.mat', kostisPath));
params_ddm_sp = array2table([params ll2bic(Gf, 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'zbias', 'bic'});

load(sprintf('%s/DDM_DC.mat', kostisPath));
params_ddm_dc = array2table([params ll2bic(Gf, 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'dcbias', 'bic'});

load(sprintf('%s/DDM_SP_DC.mat', kostisPath));
params_ddm_sp_dc = array2table([params ll2bic(Gf, 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'dcbias', 'zbias', 'bic'});

load(sprintf('%s/DDM_RP.mat', kostisPath));
params_ddm_rp = array2table([params ll2bic(Gf, 5, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'slope', 'bic'});

load(sprintf('%s/DDM_RP3.mat', kostisPath));
% So the offset is pars(7) and the slope pars(6)
params_ddm_rp2 = array2table([params ll2bic(Gf, 6, notrials)], 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'slope', 'offset', 'bic'});

params_ddm.Properties.VariableNames         = cellfun((@(x) cat(2, 'ddmK_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_z_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_dc.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_dc_', x)), params_ddm_dc.Properties.VariableNames, 'un', 0);
params_ddm_sp_dc.Properties.VariableNames   = cellfun((@(x) cat(2, 'ddmK_dcz_', x)), params_ddm_sp_dc.Properties.VariableNames, 'un', 0);
params_ddm_rp.Properties.VariableNames      = cellfun((@(x) cat(2, 'ddmK_rp_', x)), params_ddm_rp.Properties.VariableNames, 'un', 0);
params_ddm_rp2.Properties.VariableNames     = cellfun((@(x) cat(2, 'ddmK_rp2_', x)), params_ddm_rp2.Properties.VariableNames, 'un', 0);

ddm_table = cat(2, params_ddm, params_ddm_sp, params_ddm_dc, params_ddm_sp_dc, params_ddm_rp, params_ddm_rp2);

% load the main results file
results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_neutral/allindividualresults.csv');
results = results(results.session == 0, :);
results = cat(2, results, ou_table, ddm_table);
writetable(results, '/Users/urai/Data/HDDM/summary/Anke_MEG_neutral/allindividualresults_kostis.csv');

corrplot(results, {'repetition_alldata', 'ddmK_z_zbias', 'ddmK_dc_dcbias', 'ddmK_dcz_zbias', 'ddmK_dcz_dcbias'});
print(gcf, '-dpdf', '~/Data/serialHDDM/kostisData_overview_DDM.pdf');

close all;
corrplot(results, {'repetition_alldata', 'ddmK_rp_slope', 'ddmK_dc_dcbias', ...
    'ddmK_rp2_slope', 'ddmK_rp2_offset'});
print(gcf, '-dpdf', '~/Data/serialHDDM/kostisData_overview_DDM_ramp.pdf');

end