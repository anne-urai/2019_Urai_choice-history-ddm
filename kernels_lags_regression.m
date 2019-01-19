function kernels_lags_regression

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;
colors = cbrewer('qual', 'Set2', length(datasets));

set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex', ...
    'DefaultFigureWindowStyle','normal');

%% PREVIOUS STIMULUS AND RESPONSE

mdls(1).name = 'regresszlag3';
mdls(1).types{1} = {'z_prevstim', 'z_prev2stim', 'z_prev3stim'};
mdls(1).types{2} = {'z_prevresp', 'z_prev2resp', 'z_prev3resp'};

mdls(2).name = 'regressdclag3';
mdls(2).types{1} = {'v_prevstim', 'v_prev2stim', 'v_prev3stim'};
mdls(2).types{2} = {'v_prevresp', 'v_prev2resp', 'v_prev3resp'};

mdls(3).name = 'regressdczlag3';
mdls(3).types{1} = {'v_prevstim', 'v_prev2stim', 'v_prev3stim'};
mdls(3).types{2} = {'v_prevresp', 'v_prev2resp', 'v_prev3resp'};
mdls(3).types{3} = {'z_prevstim', 'z_prev2stim', 'z_prev3stim'};
mdls(3).types{4} = {'z_prevresp', 'z_prev2resp', 'z_prev3resp'};

mdls(4).name = 'regressdczlag3';
mdls(4).types{1} = {'v_prevcorrect', 'v_prev2correct', 'v_prev3correct'};
mdls(4).types{2} = {'v_preverror', 'v_prev2error', 'v_prev3error'};
mdls(4).types{3} = {'z_prevcorrect', 'z_prev2correct', 'z_prev3correct'};
mdls(4).types{4} = {'z_preverror', 'z_prev2error', 'z_prev3error'};

for m = 1:length(mdls),
    for t = 1:length(mdls(m).types),
        
        close all;
        
        subplot(4,4,1); hold on;
        
        mat = nan(length(dataset), 3);
        for d = 1:length(datasets),
            dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
            
            % recombine weights to correct and error
            for l = 1:3,
                if contains(mdls(m).types{t}{l}, 'correct')
                    dat.([mdls(m).types{t}{l} '__' mdls(m).name]) = ...
                        dat.([regexprep(mdls(m).types{t}{l}, 'correct', 'resp'), '__' mdls(m).name]) + ...
                        dat.([regexprep(mdls(m).types{t}{l}, 'correct', 'stim'), '__' mdls(m).name]);
                elseif contains(mdls(m).types{t}{l}, 'error')
                    dat.([mdls(m).types{t}{l} '__' mdls(m).name]) = ...
                        dat.([regexprep(mdls(m).types{t}{l}, 'error', 'resp'), '__' mdls(m).name]) - ...
                        dat.([regexprep(mdls(m).types{t}{l}, 'error', 'stim'), '__' mdls(m).name]);
                end
            end
            
            mat(d, :) = [nanmean(dat.([mdls(m).types{t}{1} '__' mdls(m).name])) ...
                nanmean(dat.([mdls(m).types{t}{2} '__' mdls(m).name])) ...
                nanmean(dat.([mdls(m).types{t}{3} '__' mdls(m).name]))];
            plot(1:3, mat(d, :), 'color', colors(d, :), 'linewidth', 1);
        end
        
        plot([1 3], [0 0], 'k', 'linewidth', 0.5);
        plot(1:3, nanmean(mat), 'k', 'linewidth', 1);
        [h, pval] = ttest(mat);
        if any(h>0),
            plot(find(h==1), nanmean(mat(:, (h==1))), 'ok', 'markeredgecolor', 'w', 'markerfacecolor', 'k', 'markersize', 4);
        end
        
        if contains(mdls(m).types{t}, 'stim'),
            if contains(mdls(m).types{t}, 'z'),
                ylabel('z ~ prev stim');
            elseif contains(mdls(m).types{t}, 'v'),
                ylabel('v_{bias} ~ prev stim');
            end
        elseif contains(mdls(m).types{t}, 'resp'),
            if contains(mdls(m).types{t}, 'z'),
                ylabel('z ~ prev resp');
            elseif contains(mdls(m).types{t}, 'v'),
                ylabel('v_{bias} ~ prev resp');
            end
        elseif contains(mdls(m).types{t}, 'correct'),
            if contains(mdls(m).types{t}, 'z'),
                ylabel('z ~ prev correct');
            elseif contains(mdls(m).types{t}, 'v'),
                ylabel('v_{bias} ~ prev correct');
            end
        elseif contains(mdls(m).types{t}, 'error'),
            if contains(mdls(m).types{t}, 'z'),
                ylabel('z ~ prev error');
            elseif contains(mdls(m).types{t}, 'v'),
                ylabel('v_{bias} ~ prev error');
            end
        end
        
        ylim([-0.25 0.25]);
        set(gca, 'xtick', 1:3, 'xcolor', 'k', 'ycolor', 'k');
        xlabel('Lag (past trials)');
        axis square; offsetAxes; tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_m%d_type%d.pdf', m, t));
        
    end
    
end
end