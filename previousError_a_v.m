function previousError_a_v

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global mypath datasets datasetnames colors

markers = {'d', 's', '^', 'v',  '>', '<'};
colors = cbrewer('qual', 'Set2', length(datasets));

% PREVIOUS ERROR VS CORRECT, drift rate
subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    
    dat = readtable(sprintf('%s/%s/stimcoding_prevcorrect/group_traces.csv', mypath, datasets{d}));
    
    try
        difference = dat.v_0_ - dat.v_1_;
        h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
        
    catch % plot separately for each difficulty level
        
        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1_0_$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0_0_$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_error{c})) - (dat.(driftvars_correct{c}));            
            h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues',d+(0.1*c)-0.3, 'distWidth', 0.1);
        end
    end
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k', 'ylim', [-1 1]);
ylabel({'Drift rate (v)' 'after correct - error'});
offsetAxes;

tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_driftRate.pdf'));


%% BOUNDARY SEPARATION

% PREVIOUS ERROR VS CORRECT, drift rate
close all;
sp = subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_prevcorrect/group_traces.csv', mypath, datasets{d}));
    difference = dat.a_0_ - dat.a_1_;
    h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
    
    pval = 
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'Boundary separation (a)' 'after correct - error'});
offsetAxes;
%sp.Position(1) = sp.Position(1) + 0.5;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_boundarySeparation.pdf'));


end
