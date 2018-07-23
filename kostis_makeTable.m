function kostis_makeTable

close all; clc;
global mypath datasets datasetnames
kostisPath = '~/Data/HDDM/Anke_MEG_neutral/materials_K/models';

%% 1. first take the O-U values
load(sprintf('%s/cartesius_OU/DDM.mat', kostisPath));
params_ddm = array2table(params, 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp'});

load(sprintf('%s/cartesius_OU/DDM_SP.mat', kostisPath));
params_ddm_sp = array2table(params, 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'spbias'});

load(sprintf('%s/cartesius_OU/DDM_DC.mat', kostisPath));
params_ddm_input = array2table(params, 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'inputbias'});

load(sprintf('%s/cartesius_OU/DDM_AP.mat', kostisPath));
params_ddm_lambda = array2table(params, 'variablenames', {'threshold', 'scale', 'T0', 'dv', 'bsp', 'lambdabias'});

% put into one big table
params_ddm.Properties.VariableNames = cellfun((@(x) cat(2, 'ouK_vanilla_', x)), params_ddm.Properties.VariableNames, 'un', 0);
params_ddm_sp.Properties.VariableNames = cellfun((@(x) cat(2, 'ouK_sp_', x)), params_ddm_sp.Properties.VariableNames, 'un', 0);
params_ddm_input.Properties.VariableNames = cellfun((@(x) cat(2, 'ouK_input_', x)), params_ddm_input.Properties.VariableNames, 'un', 0);
params_ddm_lambda.Properties.VariableNames = cellfun((@(x) cat(2, 'ouK_lambda_', x)), params_ddm_lambda.Properties.VariableNames, 'un', 0);

ou_table = cat(2, params_ddm, params_ddm_sp, params_ddm_input, params_ddm_lambda);

% add repetition, computed from all of Anke's trials
load(sprintf('%s/history.mat', kostisPath));
kk(6, :) = []; % remove this subject
ou_table.repetition_alldata = kk(:, 2);

% load the main results file
results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_neutral/allindividualresults.csv');
results = results(results.session == 0, :);
results = cat(2, results, ou_table);
writetable(results, '/Users/urai/Data/HDDM/summary/Anke_MEG_neutral/allindividualresults_ou.csv');

%% 2. CHECK THAT THE HISTORY BIASES (REPETITION PROBABILITY) KOSTIS USES ARE THE SAME AS IN INDIVIDUALRESULTS.CSV

% results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_neutral/allindividualresults.csv');


end