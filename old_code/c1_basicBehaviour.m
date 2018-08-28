
clear all; close all; clc;
addpath /Users/anne/Data/pupilUncertainty_FigShare/Code/Analysis
addpath /Users/anne/Dropbox/code/MEG/;
cd /Users/anne/Dropbox/code/MEG/Behaviour;
addpath('~/Documents/gramm/');

% determine the path to the data
subjectdata = subjectspecifics('ga');
results         = readtable(sprintf('%s/Data/CSV/resultsTable.csv', subjectdata.path));
results         = results(results.session > 0, :);
% make sure placebo goes first, and shorten labels for plotting
results.drug    = categorical(cellfun(@(x) x(1:4), results.drug, 'un', 0), {'plac', 'atom', 'done'}); % otherwise labels are too big
results.session = categorical(results.session); % effects coding will take the mean of both, rather than session = 0 (which doesnt exist in this data)

% design matrix:
% 1 = intercept
% 2 = session, effect coded
% 3 = drug, effect coded
% 4 = drug, also effect coded
% 5 = interaction effect 2 * 3

% ========================================== %
% SHOW A RANGE OF BEHAVIOURS
% SEPARATELY FOR DRUG AND SESSION
% INCLUDES STATISTICS (MIXED ANOVA)
% ========================================== %

allflds             = struct;

allflds(end+1).t    = 'Behaviour';
allflds(end).f      = {'accuracy', 'dprime', 'criterion', 'abscriterion', 'rt'};

allflds(end+1).t    = 'HDDM';
allflds(end).f      = {'hddm_v', 'hddm_a', 'hddm_z', 'hddm_dc', 'hddm_t'};

allflds(end+1).t    = 'Uncertainty';
allflds(end).f      = {{'pupil_error', 'pupil_correct'}, {'dprime_pupil_bin3', 'dprime_pupil_bin1'}, ...
    {'rt_error', 'rt_correct'} {'dprime_rt_bin3', 'dprime_rt_bin1'}};

allflds(end+1).t    = 'Binned correlation to dprime and criterion';
allflds(end).f      = {'pupil_dprime_corr', 'pupil_abscriterion_corr', 'rt_dprime_corr', 'rt_abscriterion_corr'};

allflds(end+1).t    = 'Repetition behaviour';
allflds(end).f      = {'repetition', 'criterionshift', 'handshift'};
results.repetition  = results.repetition - results.stimrepetition; % test

allflds(end+1).t    = 'Modulation of repetition behaviour';
allflds(end).f      = {'pupil_repetition_corr', 'pupil_criterionshift_corr', ...
    'rt_repetition_corr', 'rt_criterionshift_corr'};

allflds(end+1).t    = 'Modulation of repetition behaviour - correct vs error';
allflds(end).f      = {'pupil_correct_criterionshift_corr', 'pupil_error_criterionshift_corr', ...
    'rt_correct_criterionshift_corr', 'rt_error_criterionshift_corr'};

allflds(end+1).t    = 'Modulation of repetition behaviour';
allflds(end).f      = {{'repetition_pupil_bin3', 'repetition_pupil_bin1'}, ...
    {'criterionshift_pupil_bin3', 'criterionshift_pupil_bin1'}, ...
    {'repetition_rt_bin3', 'repetition_rt_bin1'}, {'criterionshift_rt_bin3', 'criterionshift_rt_bin1'}};

allflds(end+1).t    = 'Next trial other stuff';
allflds(end).f      = {'pupil_nextdprime_corr', 'pupil_nextabscriterion_corr', ...
    'rt_nextdprime_corr', 'rt_nextabscriterion_corr'};

allflds(end+1).t    = 'Fruend model';
allflds(end).f      = {'choiceW', 'stimW', 'pupil_choiceW', 'rt_choiceW'};

allflds(1)          = [];

for k = 1:length(allflds),
    flds = allflds(k).f;
    close all; clear g;
    
    for f = 1:length(flds),
        
        % difference between sessions, individual datapoints
        if iscell(flds{f}) % difference measure
            g(1,f) = gramm('x', results.session, 'y', results.(flds{f}{1}) - results.(flds{f}{2}),...
                'color', results.drug, 'group', results.subjnr);
            g(1,f).set_names('x', [], 'y', {sprintf('%s -', flds{f}{1}) sprintf('%s', flds{f}{2})}, ...
                'column', []);
            results.dv = results.(flds{f}{1}) - results.(flds{f}{2});
        else % one variable
            g(1,f) = gramm('x', results.session, 'y', results.(flds{f}),...
                'color', results.drug, 'group', results.subjnr);
            g(1,f).set_names('x', [], 'y', flds{f}, ...
                'column', []);
            results.dv = results.(flds{f});
        end
        
        % one line for each participant
        g(1,f).geom_line;
        g(1,f).facet_grid([], results.drug, 'force_ticks', false);
        g(1,f).axe_property('xtick', 1:2, 'xlim', [0.5 2.5]);
        g(1,f).set_text_options('base_size', 8, 'facet_scaling', 1, 'title_scaling', 0.8);
        g(1,f).geom_hline('style', 'k:'); % put a line at zero
        g(1,f).no_legend;
        %   g(1,f).set_color_options('map', 'brewer_dark');
        
        % add stats from a mixed effects ANOVA
        % each subject has a random intercept and an individual effect of session
        % test for group-level effect of session, drug, and interaction
        % effects coding: intercept is grand mean
        lme = fitlme(results, ...
            'dv ~ 1 + drug * session + (1 + session|subjnr)', ...
            'dummyvarcoding', 'effects');
        
        anov = anova(lme); % do statistics
        pvals = pval2stars(anov.pValue);
        txt = {sprintf('intercept %s', pvals{1}) ...
            sprintf('session %s, drug %s', pvals{2}, pvals{3}) ...
            sprintf('session X drug %s', pvals{4})};
        % display these stats in the title
        g(1,f).set_title(txt);
        
        % if there are several bins, show each of them
        if iscell(flds{f}),
            if ~isnan(str2double(flds{f}{1}(end))),
                g(2,f) = gramm('x', repmat([1 2 3], length(results.session), 1), 'y', ...
                    [results.([flds{f}{1}(1:end-1) '1']) results.([flds{f}{1}(1:end-1) '2']) results.([flds{f}{1}(1:end-1) '3'])],...
                    'color', results.drug, 'group', results.drug);
                g(2,f).axe_property('xtick', 1:3, 'xlim', [0.5 3.5]);
                % parse field names
                r = regexp(flds{f}{1}, '[a-z]*', 'match');
                g(2,f).set_names('Row', 'Session', 'y', r(1), 'x', r(2));
            else
                g(2,f) = gramm('x', repmat([1 2], length(results.session), 1), 'y', ...
                    [results.(flds{f}{1}) results.(flds{f}{2})],...
                    'color', results.drug, 'group', results.drug);
                r1 = regexp(flds{f}{1}, '[a-z]*', 'match');
                r2 = regexp(flds{f}{2}, '[a-z]*', 'match');
                g(2,f).axe_property('xtick', 1:2, 'xticklabel', [r1(2) r2(2)], 'xlim', [0.5 2.5]);
                g(2,f).set_names('Row', 'Session', 'y', r1(1), 'x', []);
            end
            g(2,f).stat_summary('type', 'fitnormalci', ...
                'geom', {'errorbar', 'line'}, 'setylim', true);
            g(2,f).facet_grid(results.session, []);
            g(2,f).no_legend;
            g(2,f).set_text_options('base_size', 8, 'facet_scaling', 1, 'big_title_scaling', 1);
            %   g(2,f).set_color_options('map', 'brewer_dark');
        end
    end
    g.draw;
    
    % also add SEM based on grouping variable
    for f = 1:length(flds),
        g(1,f).update('group', results.drug);
        g(1,f).stat_summary('type', 'fitnormalci', ...
            'geom', {'black_errorbar', 'line'});
        g(1,f).no_legend;
    end
    
    % add blank subplot to make sure subplots are not super thin
    g(3,f-1) = gramm('x', results.session, 'y', results.subjnr,...
        'color', results.drug, 'group', results.subjnr);
    g(3,f-1).axe_property('visible', 'off');
    g(3,f-1).no_legend;
    g.draw;
    
    % make the lines black
    for f = 1:length(flds),
        set([g(1,f).results.stat_summary.line_handle], 'Color','k');
    end
    
    % save this figure
    suptitle(allflds(k).t); % which type of behaviour is this?
    set(gcf, 'PaperPositionMode', 'auto'); % avoid a warning
    print(gcf, '-dpdf', sprintf('%s/Figures/behaviour_%d.pdf', subjectdata.path, k));
end

% ========================================== %
% BETWEEN-SUBJECT CORRELATIONS
%% ========================================== %

clearvars -except results subjectdata

allflds             = struct;
allflds(end+1).t    = 'Criterion shift';
allflds(end).f      = {{'criterionshift', 'pupil_criterionshift_corr'}, {'criterionshift', 'rt_criterionshift_corr'}};

allflds(end+1).t    = 'Repetition';
allflds(end).f      = {{'repetition', 'pupil_repetition_corr'}, {'repetition', 'rt_repetition_corr'}};
allflds(1)          = [];

for k = 1:length(allflds),
    clf; clear g;
    
    for f = 1:length(allflds(k).f),
        % first show the linear fit
        g(1,f) = gramm('x', results.(allflds(k).f{f}{1}), ...
            'y', results.(allflds(k).f{f}{2}), 'color', results.drug);
        
        % manual stats
        corrFunc = @(x,y) corr(x,y, 'rows', 'complete');
        [rho1, pval1] = splitapply(corrFunc, results.(allflds(k).f{f}{1})(double(results.session) == 1), ...
            results.(allflds(k).f{f}{2})(double(results.session) == 1), findgroups(results.drug(double(results.session) == 1)));
        [rho2, pval2] = splitapply(corrFunc, results.(allflds(k).f{f}{1})(double(results.session) == 2), ...
            results.(allflds(k).f{f}{2})(double(results.session) == 2), findgroups(results.drug(double(results.session) == 2)));
        
        [~, n] = findgroups(results.drug);
        txt = {sprintf('%s S1 r = %.3f p = %.3f, S2 r = %.3f p = %.3f', char(n(1)), rho1(1), pval1(1), rho2(1), pval2(1)) ...
            sprintf('%s S1 r = %.3f p = %.3f, S2 r = %.3f p = %.3f', char(n(2)), rho1(2), pval1(2), rho2(2), pval2(2)) ...
            sprintf('%s S1 r = %.3f p = %.3f, S2 r = %.3f p = %.3f', char(n(3)), rho1(3), pval1(3), rho2(3), pval2(3))};
        g(1,f).set_title(txt);
        g(1,f).set_text_options('base_size', 8, 'facet_scaling', 1, 'title_scaling', 1);
        g(1,f).stat_glm(); % show regression line
        %g(1,f).set_color_options('map', [0 0 0]); % black
        g(1,f).no_legend();
    end
    g.draw;
    
    for f = 1:length(allflds(k).f),
        % one point for each participant
        g(1,f).update('color', results.drug);
        g(1,f).geom_point;
        g(1,f).facet_grid([], results.session, 'force_ticks', false);
        g(1,f).set_text_options('base_size', 8, 'facet_scaling', 1, 'title_scaling', 0.8);
        g(1,f).set_color_options('map', 'lch');
        g(1,f).set_names('x', allflds(k).f{f}{1}, 'y', allflds(k).f{f}{2}, 'column', 'Session');
        g(1,f).no_legend();
        g(1,f).axe_property('PlotBoxAspectRatio', [1 1 1]); % square
    end
    g.draw;
    
    suptitle(allflds(k).t); % which type of behaviour is this?
    print(gcf, '-dpdf', sprintf('%s/Figures/behaviour_correlations_%d.pdf', subjectdata.path, k));
end


