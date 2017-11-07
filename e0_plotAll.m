
%% determine how the figures will look
clear all; clc; close all;
set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');

global datasets datasetnames mypath

usr = getenv('USER');
switch usr
    case 'anne'
        
        mypath = '~/Data/HDDM';
    case 'aeurai'
        mypath  = '/nfs/aeurai/HDDM';
end

% neutral vs biased plots
datasets = {'RT_RDK', 'JW_yesno', 'NatComm', 'MEG', ...
    'JW_PNAS', 'JW_fMRI', ...
    'Anke_2afc_sequential', ...
    'MEG_MEGsessions',  'Bharath_fMRI', 'Anke_MEG', 'Anke_merged'};
datasets = {'Murphy', 'JW_yesno', 'NatComm', 'MEG', 'JW_PNAS', ...
    'Anke_2afc_sequential', 'Anke_MEG', 'Bharath_fMRI'};

datasetnames = { {'2AFC motion RT'},  {'Yes/no tone RT'}, ...
    {'2IFC-1'}, {'2IFC-2'}, ...
    {'Yes/no contrast RT'},  ...
    {'2AFC-1, Braun et al. 2017'}, ...
   {'2AFC-2, Anke MEG'}, {'2AFC-3, Bharath fMRI'}};

%% datasets = datasets(7);
% datasetnames = datasetnames(7);
% go to code
try
    cd('/Users/anne/Drive/Dropbox/code/RT_RDK');
end

%% start the actual plots

disp('starting');
% f0_schematic_DDM_bias; % figure 3a

e6_serialBias_SfN_modelFree_CRF_PPC

% sv_comparison;
e2_serialBias_SfN_SanityChecks; % correlate dprime with drift rate
e8_serialBias_SfN_PPC; % figure 2, show that all models fit OK

e1_serialBias_SfN_DIC; % figure 3b & c
%e1_serialBias_SfN_BIC;

% % show the fits separately for dc and z
% alldat = e1b_serialBias_SfN_ModelFreeCorrelation_independentFits; % figure 4
% forestPlot(alldat);
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_indep.pdf'));

% main figure: correlations from the jointly fit model
close all;
for Gsq = [0],
    for sz = [0],
        
        if Gsq == 0 && sz == 1, continue; end % hierarchical sampling with sz takes forever
        
        alldat = e1b_serialBias_SfN_ModelFreeCorrelation_grey(Gsq, sz); % figure 4
        forestPlot(alldat);
        
        switch Gsq
            case 1
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_sz%d_Gsq.pdf', sz));
            case 0
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_sz%d_HDDM.pdf', sz));
        end
    end
end


% ======================= %
% REGRESSION MODELS
% ======================= %

