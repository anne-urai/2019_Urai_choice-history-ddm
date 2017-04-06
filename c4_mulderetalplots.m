
clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 1:length(datasets),
    
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));
    
    % quick plot of mulder et al.
    cnt = 1; figure;
    sessions = [unique(results.session)'];
    sessions(sessions == 0) = []; % dont use the stuff computed over all sessions
    
    for s = sessions,
        
        subplot(5,5,[cnt:cnt+1]); cnt = cnt + 3;
        hold on;
        
        % find only repeaters for now
        idx = find(results.session == s);
        
        % first, correct
        dat = [results.rt_invalid_fast_correct(idx) ...
            results.rt_invalid_slow_correct(idx) results.rt_valid_slow_correct(idx)...
            results.rt_valid_fast_correct(idx)];
        %plot(1:2, dat, '.-g');
        errorbar(1:4, nanmean(dat), nanstd(dat) ./ sqrt(size(dat, 1)), 'g');
        
        %pval = signrank(dat(:, 1), dat(:, 4));
        %mysigstar(gca, [1 4], max(get(gca, 'ylim')), pval);
        
        % now errors in red
        dat = [results.rt_invalid_fast_error(idx) ...
            results.rt_invalid_slow_error(idx) results.rt_valid_slow_error(idx)...
            results.rt_valid_fast_error(idx)];
        
        %plot(3:4, dat, '.-r');
        errorbar(gca, 5:8, nanmean(dat), nanstd(dat) ./ sqrt(size(dat, 1)), 'r');
        %pval = signrank(dat(:, 1), dat(:, 5));
        %mysigstar(gca, [5 8], max(get(gca, 'ylim')), pval);
        ylabel('RT (s)');
        
        set(gca, 'xtick', 1:8, 'xticklabel', ...
            repmat({'invalid fast', 'invalid slow', 'valid slow', 'valid fast'}, 1, 2), ...
            'xlim', [0.5 8], 'xticklabelrotation', -30);
        axisNotSoTight;
        title(sprintf('Session %d', s));
        if s < max(results.session), set(gca, 'xticklabel', []); end
        
        % second y axis
        subplot(5,5,cnt); cnt = cnt + 2;
        
        dat = [results.accuracy_invalid_fast(idx) ...
            results.accuracy_invalid_slow(idx) results.accuracy_valid_slow(idx)...
            results.accuracy_valid_fast(idx)];
        %plot(1:2, dat, '.-g');
        errorbar(gca, 1:4, nanmean(dat), nanstd(dat) ./ sqrt(size(dat, 1)), 'k');
        ylabel('Accuracy')
        %pval = signrank(dat(:, 1), dat(:, 2));
        %mysigstar(gca, [1 2], max(get(gca, 'ylim')), pval);
        box off; axisNotSoTight;
        
        set(gca, 'xtick', 1:4, 'xticklabel', ...
            {'invalid fast', 'invalid slow', 'valid slow', 'valid fast'}, ...
            'xlim', [0.5 4], 'xticklabelrotation', -30);
        if s < max(results.session), set(gca, 'xticklabel', []); end
        
    end
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/mulderetalplot_prevresp.pdf', datasets{d}));
    
end

% pl