% get a huge list with values for each participant
% can then work with this dataframe
clear all; close all; clc;
addpath /Users/anne/Data/pupilUncertainty_FigShare/Code/Analysis
addpath /Users/anne/Dropbox/code/MEG/;
cd /Users/anne/Dropbox/code/MEG/Behaviour;
subjectdata    = subjectspecifics('ga');

% compute a bunch of basic things from Matlab
results        = A5a_behaviouralMetrics();
close all;

% add in the Fruend weights from modelfits
results        = A5b_FruendModel(results);

% add HDDM fits!
params = {'a', 'v', 't', 'z', 'dc'};
% boundary separation, drift rate, non-decisiontime, starting point, drift criterion
results.hddm_a  = nan(size(results.subjnr));
results.hddm_v  = nan(size(results.subjnr));
results.hddm_t  = nan(size(results.subjnr));
results.hddm_z  = nan(size(results.subjnr));
results.hddm_dc = nan(size(results.subjnr));

for md = 1,
    hddm = readtable(sprintf('%s/Data/HDDM/2ifc_MEGdata_1/diagnostics/results%d.csv', ...
        subjectdata.path, md));
    for p = 1:length(params),
        for sj = subjectdata.all,
            thissubjectdata = subjectspecifics(sj);
            for s = 1:length(thissubjectdata.session),
                results.(sprintf('hddm_%s', params{p})) ...
                    (results.subjnr == sj & results.session == s) = ...
                    hddm{strcmp(hddm.Var1, sprintf('%s_subj(%s.%d).%d', ...
                    params{p}, thissubjectdata.drug, s, sj)), 2};
            end
        end
    end
end
writetable(results, sprintf('%s/Data/CSV/resultsTable.csv', subjectdata.path));

if 0,
% ========================================== %
% SANITY CHECKS
% ========================================== %

flds2 = {'pupil', 'rt'};
flds1 = {'repetitioncrit', 'criterionshift', 'repetition', 'stimrepetition'};
cnt = 0;
for f1 = 1:2,
    for f2 = 1:2,
        cnt = cnt + 1;
        subplot(2,2,cnt);
        scatter([results.(sprintf('%s_%s_bin%d', flds1{f1}, flds2{f2}, 3))(results.session == 0) - ...
            results.(sprintf('%s_%s_bin%d', flds1{f1}, flds2{f2}, 1))(results.session == 0)], ...
            results.(sprintf('%s_%s_corr', flds2{f2}, flds1{f1}))(results.session == 0), ...
            20, findgroups(results.drug(results.session == 0)), 'filled');
        grid on;
        l = lsline; l.Color = 'k';
        axis tight;  axis square; axisEqual;
        r = refline(1); r.Color = [0.5 0.5 0.5];
        title(sprintf('%s, %s', flds2{f2}, flds1{f1}));
    end
end
suplabel('High - Low', 'x'); suplabel('Correlation', 'y');

% ========================================== %
% plot stuff bewteen sessions
% ========================================== %

set(groot, 'DefaultAxesFontSize', 7);
figure;
flds = {'hddm_a', 'hddm_v', 'hddm_t', 'hddm_z', 'hddm_dc', ...
    'dprime', 'accuracy', 'criterion', 'criterionshift'};
colors = viridis(3);
nsubpl = ceil(sqrt(length(flds)));
for f = 1:length(flds),
    subplot(nsubpl,nsubpl,f);
    
    h = gscatter(results.(flds{f})(results.session == 1), ...
        results.(flds{f})(results.session == 2),  ...
        (results.drug(results.session == 1)), ...
        colors(1:3, :), '.', 10, 0);
    
    ylabel(flds{f}, 'interpreter', 'none');
    axis square; box off; axis tight;
    ylims = get(gca, 'ylim'); xlims = get(gca, 'xlim');
    uselims = [min([ylims(1) xlims(1)]) max([ylims(2) xlims(2)])];
    ylim(uselims); xlim(uselims);
    
    [rho, pval] = corr(results.(flds{f})(results.session == 1), ...
        results.(flds{f})(results.session == 2), ...
        'type', 'spearman', 'rows', 'pairwise');
    if pval < (0.01), l = lsline; end
    xlabel(sprintf('r=%.3f p=%.3f', rho, pval), 'fontweight', 'normal');
    grid on;
    
end

try
    s = subplot(5,5,f+1);  spos  = s.Position;
    [d, g] = findgroups(results.drug(results.session == 1));
    l = legend(h, g);
    l.Position = spos; l.Box = 'off'; axis off;
end
print(gcf, '-dpdf', sprintf('%s/Figures/dataFrameOverview.pdf', subjectdata.path));

% ========================================== %
% compare fruend weights with other metrics
% ========================================== %

fruendVars = {'choiceW', 'stimW', 'pupil_choiceW', 'pupil_stimW', 'rt_choiceW', 'rt_stimW'};
matlVars   = {'pupil_repetition',  'pupil_repetitioncrit','pupil_criterionshift',  ...
    'rt_repetition', 'rt_repetitioncrit', 'rt_criterionshift' };
corrplot(results(results.session == 0, :), fruendVars, matlVars);
end
