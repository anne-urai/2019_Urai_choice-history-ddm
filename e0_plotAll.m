
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
datasets = {'Murphy', 'JW_yesno', 'NatComm', 'MEG', 'Anke_2afc_sequential', ...
   'MEG_MEGsessions',  'Bharath_fMRI', 'Anke_MEG', 'Anke_merged'};

datasetnames = { {'2AFC RT'},  {'Yes/no RT'}, ...
    {'2IFC-1'}, {'2IFC-2'}, ...
    {'2AFC-1, Braun et al. 2017'}, ...
    {'2IFC-2, MEG sessions'}, ...
    {'2AFC-2'}, {'2AFC-1, Anke MEG'}, {'2AFC-1'}};

% go to code
try
cd('/Users/anne/Drive/Dropbox/code/RT_RDK');
end

%% start the actual plots

disp('starting');
% f0_schematic_DDM_bias; % figure 3a

% e2_serialBias_SfN_SanityChecks; % correlate dprime with drift rate
%e8_serialBias_SfN_PPC; % figure 2, show that all models fit OK

%e1_serialBias_SfN_DIC; % figure 3b & c

alldat = e1b_serialBias_SfN_ModelFreeCorrelation_independentFits; % figure 4
forestPlot(alldat);
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_indep.pdf'));

close all;
alldat = e1b_serialBias_SfN_ModelFreeCorrelation_grey; % figure 4
forestPlot(alldat);
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot.pdf'));
print(gcf, '-depsc', sprintf('~/Data/serialHDDM/forestplot.eps'));

%alldat = e1b_serialBias_SfN_ModelFreeCorrelation_chiSquare; % figure 4
%forestPlot(alldat);
%print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_chiSquare.pdf'));

%pharmaPosteriors;
%e9_prevCorrect; % error vs correct
%e3_serialBias_SfN_Posteriors; % figure 4c

% to add (from regression models)
% - RT modulation
% - kernels (lags)

%
