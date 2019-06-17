function post_error_slowing

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global mypath datasets datasetnames 
colors = cbrewer('qual', 'Set2', length(datasets));

% ========================================= %
% POST-ERROR SLOWING
% ========================================= %

close all;
sp = subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    
    % COMPUTE POST ERROR SLOWING
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    PES = results.posterrorslowing;
    h = scatter(d*ones(size(PES)), PES, 3, colors(d, :), 'jitter', 'on', 'jitteramount', 0.05);
    plot([d-0.2 d+0.2], [nanmean(PES) nanmean(PES)], 'k-');
    pval = permtest(PES);
    mysigstar(gca, d, -0.13, pval);
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
     
end

ylim([-0.15 0.25]);
set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'Post-error slowing' 'RT_{e+1} - RT_{c+1} (s)'});
offsetAxes;
xlim([ -0.25 7]);
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PostErrorSlowing_traditional.pdf'));

% ========================================= %
% PREVIOUS ERROR VS CORRECT, drift rate
% ========================================= %

close all;
subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    
    dat = readtable(sprintf('%s/%s/stimcoding_prevcorrect/group_traces.csv', mypath, datasets{d}));
    
    switch datasets{d}
    case 'NatComm'
        
        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1_0_$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0_0_$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_error{c})) - (dat.(driftvars_correct{c}));            
           % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues',d+(0.1*c)-0.3, 'distWidth', 0.1, 'histOpt', 1.1);
		    violinPlot_distribution(d+(0.15*c)-0.4, difference, colors(d, :), 25);
			
        end

    case 'Anke_MEG_transition'

        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1_$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0_$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_error{c})) - (dat.(driftvars_correct{c}));            
           % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues',d+(0.1*c)-0.3, 'distWidth', 0.1, 'histOpt', 1.1);
            violinPlot_distribution(d+(0.15*c)-0.4, difference, colors(d, :), 25);
        end
        
    otherwise
        difference = dat.v_0_ - dat.v_1_;
       % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
        violinPlot_distribution(d, difference, colors(d, :));

    end
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k', 'ylim', [-1 1]);
ylabel({'Drift rate (v)' 'after error - correct'});
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_driftRate.pdf'));

% ========================================= %
% BOUNDARY SEPARATION
% ========================================= %

close all;
sp = subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_prevcorrect/group_traces.csv', mypath, datasets{d}));
    difference = dat.a_0_ - dat.a_1_;
    % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
    violinPlot_distribution(d, difference, colors(d, :));
	
end
ylim([-0.6 0.6]);

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'Boundary separation (a)' 'after error - correct'});
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_boundarySeparation.pdf'));






% ================================================================================== %
% NOW DO THIS AGAIN - BUT FOR MODEL WITH ALSO CHOICE HISTORY BIAS!
% ================================================================================== %

close all;
subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_PES/group_traces.csv', mypath, datasets{d}));
    
    switch datasets{d}
    case 'NatComm'
        
        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1_0_$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0_0_$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_error{c})) - (dat.(driftvars_correct{c}));            
           % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues',d+(0.1*c)-0.3, 'distWidth', 0.1, 'histOpt', 1.1);
            violinPlot_distribution(d+(0.15*c)-0.4, difference, colors(d, :), 25);
            
        end

    case 'Anke_MEG_transition'

        vars    = dat.Properties.VariableNames';
        driftvars_correct   = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_1_$')));
        driftvars_error     = vars(~cellfun(@isempty, regexp(vars, 'v_\S+_0_$')));
        
        for c = 1:length(driftvars_correct),
            difference = (dat.(driftvars_error{c})) - (dat.(driftvars_correct{c}));            
           % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues',d+(0.1*c)-0.3, 'distWidth', 0.1, 'histOpt', 1.1);
            violinPlot_distribution(d+(0.15*c)-0.4, difference, colors(d, :), 25);
        end
        
    otherwise
        difference = dat.v_0_ - dat.v_1_;
       % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
        violinPlot_distribution(d, difference, colors(d, :));

    end
end

set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k', 'ylim', [-1 1]);
ylabel({'Drift rate (v)' 'after error - correct'});
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_choiceHist_driftRate.pdf'));

% ========================================= %
% BOUNDARY SEPARATION
% ========================================= %

close all;
sp = subplot(3,3,1);
hold on;
% plot identity line
plot([1 6], [0 0], 'k-', 'linewidth', 0.5);

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/%s/stimcoding_dc_z_PES/group_traces.csv', mypath, datasets{d}));
    difference = dat.a_0_ - dat.a_1_;
    % h = violinPlot(difference, 'color', colors(d, :), 'showMM', 6, 'xValues', d);
    violinPlot_distribution(d, difference, colors(d, :));
    
end

ylim([-0.6 0.6]);
set(gca, 'xtick', 1:length(datasets), 'xticklabel', legtext, 'xticklabelrotation', -30, 'xcolor', 'k');
ylabel({'Boundary separation (a)' 'after error - correct'});
offsetAxes;
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevErrorCorrect_PES_choiceHist_boundarySeparation.pdf'));



end
