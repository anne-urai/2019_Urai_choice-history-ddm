
%% determine how the figures will look
clear all; clc; close all;
set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');

global datasets datasetnames mypath colors

usr = getenv('USER');
switch usr
    case 'anne'
        mypath = '~/Data/HDDM';
    case 'aeurai'
        mypath  = '/nfs/aeurai/HDDM';
end

% neutral vs biased plots
datasets = {'Murphy', 'JW_yesno', 'JW_PNAS', 'NatComm', 'MEG'}; 
datasetnames = { {'Left/right motion RT'},  {'Yes/no tone RT'}, {'Yes/no contrast RT'}, ...
    {'2IFC-1'}, {'2IFC-2'}, {'Up/down motion RT'}, ...
    {'2AFC-1, Braun et al. 2017'}, ...
   {'2AFC-2, Anke MEG'}, {'2AFC-3, Bharath fMRI'}};

% go to code
try
    cd('/Users/anne/Drive/Dropbox/code/RT_RDK');
end

% from Thomas, green blue grey
colors = [141 165 8;  8 141 165; 150 150 150] ./ 256;
colors = [51,160,44; 31,120,180] ./ 256;
colors = [178,223,138; 166,206,227] ./ 256; % lighter
colors = [77,175,74; 55,126,184] ./ 256; % green blue

%% PREPARING DATA
if 0,
	b2_HDDM_readIntoMatlab(datasets);
	b2b_Gsq_readIntoMatlab(datasets);
	b3_makeDataframe(datasets);
end

disp('starting');

% ======================= %
% SANITY CHECKS/ MODEL FITS
% ======================= %

 sv_comparison;
 e3_serialBias_SfN_repetitionRange
 e2_serialBias_SfN_SanityChecks; % correlate dprime with drift rate
 e1_serialBias_SfN_DIC; % figure 3b & c
 e8_serialBias_SfN_PPC; % figure 2, show that all models fit OK
 e1_serialBias_SfN_BIC;

% % show the fits separately for dc and z
% alldat = e1b_serialBias_SfN_ModelFreeCorrelation_independentFits; % figure 4
% forestPlot(alldat);
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_indep.pdf'));

% ======================= %
% MODEL FREE CONFIRMATION
% ======================= %

e6_serialBias_SfN_modelFree_CRF_PPC

% ======================= %
% CORRELATIONS WITH P(REPEAT)
% ======================= %

close all;
for sz = [0 1],
	for Gsq = [0 1],
        
         % if Gsq == 0 && sz == 1, continue; end % hierarchical sampling with sz takes forever
        
        alldat = e1b_serialBias_SfN_ModelFreeCorrelation_grey(Gsq, sz); % figure 4
        forestPlot(alldat);
        
        switch Gsq
            case 1
                filename = sprintf('~/Data/serialHDDM/forestplot_sz%d_Gsq.pdf', sz);
            case 0
                filename = sprintf('~/Data/serialHDDM/forestplot_sz%d_HDDM.pdf', sz);
        end
		print(gcf, '-dpdf', filename);
		disp(filename);
    end
end


% ======================= %
% PREVCORRECT
% ======================= %

% alldat = e1b_serialBias_SfN_ModelFreeCorrelation_prevCorrect;
% forestPlot(alldat);
% print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_HDDM_prevcorrect.pdf'));

% ======================= %
% REGRESSION MODELS
% ======================= %

% 1. compare dic between models with just history and with neural data
%e11_serialBias_MEGregression_DIC
%e11_serialBias_MEGregression_posteriors

% ========================= %
% MEG PHARMA
% ========================= %

alldat = e1b_serialBias_SfN_ModelFreeCorrelation_MEGpharma(); % figure 4
forestPlot(fliplr(alldat));
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/forestplot_pharma.pdf'));


% ======================= %
% SCHEMATIC/HYPOTHESES
% ======================= %

 f0_schematic_DDM_bias; % figure 3a

