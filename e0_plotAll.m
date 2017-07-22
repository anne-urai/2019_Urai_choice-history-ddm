
%% determine how the figures will look
set(groot, 'defaultaxesfontsize', 6.5, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1.2, ...
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
        mypath  = '~/Data/HDDM';
end

% neutral vs biased plots
datasets = {'RT_RDK', 'NatComm', 'MEG', 'Anke_2afc_alternating', 'Anke_2afc_neutral', 'Anke_2afc_repetitive'};
datasetnames = { {'2AFC, RT', 'n = 22'}, ...
    {'2IFC, Urai et al. 2017', 'n = 27'}, {'2IFC, replication', 'n = 61'}, ...
    {'2AFC, Braun et al. 2017', 'Alternating'}, ...
    {'2AFC, Braun et al. 2017', 'Neutral'}, ...
    {'2AFC, Braun et al. 2017', 'Repetitive'}};

%% start the actual plots

disp('starting');
e8_serialBias_SfN_PPC; % figure 2, show that all models fit OK

assert(1==0);
f0_schematic_DDM_bias; % figure 3a
e1_serialBias_SfN_DIC; % figure 3b & c
e1b_serialBias_SfN_ModelFreeCorrelation; % figure 4
e3_serialBias_SfN_Posteriors; % figure 4c
e9_prevCorrect; % error vs correct

% to add (from regression models)
% - RT modulation
% - kernels (lags)

%
