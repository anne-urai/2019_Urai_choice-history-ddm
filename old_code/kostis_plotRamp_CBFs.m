function kostis_plotRamp_CBFs
% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
global mypath datasets datasetnames colors

%% TODO: IN FITS, SAVE ALL RTS AND CHOICES TO COMPUTE CONDITIONAL BIAS FUNCTIONS