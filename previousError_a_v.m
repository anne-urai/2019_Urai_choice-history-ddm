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

markers = {'o', '^', 'v', 's', 'd', '+'};
colors = cbrewer('qual', 'Set2', length(datasets));

% PREVIOUS ERROR VS CORRECT, drift rate
subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    try
        difference = dat.v_1__stimcodingprevcorrect - dat.v_0__stimcodingprevcorrect;
        
        h = ploterr(d, nanmean(difference), [], nanstd(difference) ./ sqrt(length(find(~isnan(difference)))), ...
            markers{d}, 'abshhxy', 0);
        set([h(1) h(2)], 'color', colors(d, :));
        set(h(1), 'markerfacecolor', 'w', 'markersize', 4);
    catch % plot separately for each difficulty level
        
        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1__stimcodingprevcorrect$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0__stimcodingprevcorrect$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_correct{c})) - (dat.(driftvars_error{c}));
            h = ploterr(d+(0.1*c)-0.3, nanmean(difference), [], nanstd(difference) ./ sqrt(length(find(~isnan(difference)))), ...
                markers{d}, 'abshhxy', 0);
            set([h(1) h(2)], 'color', colors(d, :));
            set(h(1), 'markerfacecolor', 'w', 'markersize', 2+0.5*c);
        end
    end
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30);
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
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    difference = dat.a_1__stimcodingprevcorrect - dat.a_0__stimcodingprevcorrect;
    
    h = ploterr(d, nanmean(difference), [], nanstd(difference) ./ sqrt(length(find(~isnan(difference)))), ...
        markers{d}, 'abshhxy', 0);
    set([h(1) h(2)], 'color', colors(d, :));
    set(h(1), 'markerfacecolor', 'w', 'markersize', 4);
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30);
ylabel({'Boundary separation (a)' 'after correct - error'});
offsetAxes;
%sp.Position(1) = sp.Position(1) + 0.5;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_boundarySeparation.pdf'));


end
