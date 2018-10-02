% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

%% ========================================== %
% determine how the figures will look
% ========================================== %

clear all; clc; close all;
set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex', ...
    'DefaultFigureWindowStyle','normal');

global datasets datasetnames mypath colors

usr = getenv('USER');
switch usr
    case {'anne', 'urai'}
        mypath = '~/Data/HDDM';
    case 'aeurai'
        mypath  = '/nfs/aeurai/HDDM';
end

% neutral vs biased plotsC
% reorder the datasets: first show the ?standard ones? ?visual motion 2AFC with
% feedback, then 2IFC, and only then yes-no / no feedback.
datasets = {'Murphy', 'Anke_MEG_transition', 'NatComm', 'MEG', 'JW_PNAS', 'JW_yesno', 'Anke_MEG_transition_no81'};

datasetnames = {{'Visual motion' '2AFC (RT)'},   {'Visual motion' '2AFC (FD)'},...
    {'Visual motion' '2IFC (FD) #1'}, {'Visual motion' '2IFC (FD) #2'}, ...
    {'Visual contrast' 'yes/no (RT)'}, {'Auditory' 'yes/no (RT)'}, ...
     {'Visual motion' '2AFC (FD)'}};

% go to code
try
    addpath(genpath('~/code/Tools'));
    cd('/Users/urai/Documents/code/serialDDM');
end

% from Thomas, green blue grey
colors = [77,175,74; 55,126,184] ./ 256; % green blue

%% ========================================== %
% PREPARING DATA
% This will generate the allindividualresults.csv files
% ========================================== %

%datasetnames = datasetnames([2 7]);
% datasets = datasets([2 7]);
if 0,
    read_into_Matlab(datasets);
    read_into_Matlab_gSquare(datasets);
    make_dataframe(datasets);
    % rename_PPC_files(datasets);
end

disp('starting');
return;

% ========================================== %
% Figure 1. SCHEMATIC/HYPOTHESES
% ========================================== %

schematic;

% ========================================== %
% FIGURE 2
% ========================================== %

repetition_range;
strategy_plot;

% ========================================== %
% FIGURE 3
% ========================================== %

barplots_DIC; %3a
 
close all; subplot(3,3,1);
conditional_bias_functions_collapsed(4, 3, 'rt', 0, 0);
tightfig; print(gcf, '-dpdf', '~/Data/serialHDDM/conditionalBiasFunctions_collapsed.pdf'); % 3b

% barplots
close all; conditional_bias_functions_collapsed_summary;
print(gcf, '-dpdf', '~/Data/serialHDDM/conditionalBiasFunctions_barplots.pdf'); % 3c

% ========================================== %
% FIGURE 4
% ========================================== %

alldat = individual_correlation_main(0, 0); % figure 4
forestPlot(alldat);
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_st%d_HDDM.pdf', 0));

% ========================================== %
% FIGURE 5
% ========================================== %

alldat = individual_correlation_prevcorrect;
% separate plots for correct and error
forestPlot(alldat(1:2:end));
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_HDDM_prevcorrect.pdf'));
forestPlot(alldat(2:2:end));
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_HDDM_preverror.pdf'));

% compare the correlation coefficients for figure 5d
compare_correlations_correct_error(alldat);

% DIC comparison
barplots_DIC_previousresponse_outcome;

% ========================================== %
% FIGURE 6
% ========================================== %

% grab the results from Kostis' fits
% kostis_driftRate;
kostis_makeTable_v2;

% a. motion energy filtering: this has to run before, 
% on the UKE cluster to grab Anke's motionenergy coordinates
motionEnergy_filterDots;
motionEnergy_check;
motionEnergy_kernels;

% a. normal DDM
kostis_plotDDM_BIC;
kostis_plotDDM_correlation;

% b. ramping vs. static drift bias
kostis_plotRamp_correlation;
kostis_plotRamp_BIC;

% c. DDM with collapsing bounds
kostis_plotDDMCol_BIC;
kostis_plotDDMCol_correlation;


% PART 2: O-U FITS
kostis_plotOU_BIC;
kostis_plotOU_correlation;

% PART 2: O-U FITS
kostis_plotOUD_BIC;
kostis_plotOUD_correlation;

% ========================================== %
% SUPPLEMENTARY FIGURE 1
% ========================================== %

% see graphicalModels.manualGraphical.py
% run in Python: plot_HDDM_priors.py

% ========================================== %
% SUPPLEMENTARY FIGURE 2
% ========================================== %

dprime_driftrate_correlation;
posterior_predictive_checks;
history_kernels;
strategy_plot_2-7;

% ========================================== %
% SUPPLEMENTARY FIGURE 3
% ========================================== %

plot_posteriors;

% ========================================== %
% SUPPLEMENTARY FIGURE 4
% ========================================== %

% a. G-square fit
alldat = individual_correlation_main(1, 0);
forestPlot(alldat);
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_st%d_Gsq.pdf', 0));
barplots_BIC;

% b. fit with between-trial variability in non-decision time
alldat = individual_correlation_main(0, 1); %
forestPlot(alldat);
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_st%d_HDDM.pdf', 1));

% ========================================== %
% SUPPLEMENTARY FIGURE 5
% ========================================== %

alldat = individual_correlation_pharma();
forestPlot(fliplr(alldat));
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_pharma.pdf'));

% ========================================== %
% SUPPLEMENTARY FIGURE 6
% ========================================== %

post_error_slowing;

% ========================================== %
% SUPPLEMENTARY FIGURE 7
% ========================================== %



% ========================================== %
% SUPPLEMENTARY FIGURE 8
% ========================================== %

multiplicative_vbias_psychfuncs_ppc;
multiplicative_vbias_DIC;

% ========================================== %
% SUPPLEMENTARY FIGURE 9
% see JW's code in simulations/ folder
% ========================================== %

