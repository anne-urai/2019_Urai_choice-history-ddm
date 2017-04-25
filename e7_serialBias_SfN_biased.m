% ======================================================= %
% ANKE'S DATA - which parameter changes with transprob?
% ======================================================= %

% test that z indeed means more repetition in stimcoding

close all; clear; clc;

traces = readtable(sprintf('~/Data/Anke_2afc_sequential/HDDM/regress_dc_z_prevresp_prevstim/group_traces.csv'));
colors = flipud(linspecer(3));

% correlate higher z with more repetition
z_link_func = @(x) 1./(1+exp(x))';
trvars = traces.Properties.VariableNames;
zcells = find(~cellfun(@isempty, strfind(trvars, 'z_')));
for z = zcells,
    traces.(trvars{z}) = z_link_func(traces.(trvars{z}))';
end

% recode some stuff
params = {'v', 'z'};
for p = 1:length(params),
    tp = {'0_20000000000000001', '0_5', '0_80000000000000004'};
    for t = 1:length(tp),
        traces.(sprintf('%s_prevcorrect_C_transitionprob__%s_', params{p}, tp{t})) = ...
            traces.(sprintf('%s_prevresp_C_transitionprob__%s_', params{p}, tp{t})) + ...
            traces.(sprintf('%s_prevstim_C_transitionprob__%s_', params{p}, tp{t}));
        traces.(sprintf('%s_preverror_C_transitionprob__%s_', params{p}, tp{t})) = ...
            traces.(sprintf('%s_prevresp_C_transitionprob__%s_', params{p}, tp{t})) - ...
            traces.(sprintf('%s_prevstim_C_transitionprob__%s_', params{p}, tp{t}));
    end
end

params = {'v_prevresp', 'v_prevstim', 'z_prevresp', 'z_prevstim', ...
    'v_prevcorrect', 'v_preverror', 'z_prevcorrect', 'z_preverror' };
for p = 1:length(params),
    subplot(4,4,p); hold on;
    histogram(traces.(sprintf('%s_C_transitionprob__0_20000000000000001_', params{p})), ...
        'displaystyle', 'stairs', 'edgecolor', colors(1, :), 'normalization', 'pdf');
    histogram(traces.(sprintf('%s_C_transitionprob__0_5_', params{p})), ...
        'displaystyle', 'stairs', 'edgecolor', colors(2, :), 'normalization', 'pdf');
    histogram(traces.(sprintf('%s_C_transitionprob__0_80000000000000004_', params{p})), ...
        'displaystyle', 'stairs', 'edgecolor', colors(3, :), 'normalization', 'pdf');
    xlabel(regexprep(params{p}, '\_', ' ~ '));
    
end

l             = legend({'Repetitive', 'Neutral', 'Alternating'});
l.Position(2) = l.Position(2) - 0.25;
l.Box         = 'off';

print(gcf, '-dpdf', '~/Data/serialHDDM/fig7_Anke_conditions.pdf');

% ======================================================= %
% ANKE'S DATA - correlate HDDM with p(repeat)
% ======================================================= %

close all; clear; clc;
dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
    'Anke_2afc_sequential'));

% correlate higher z with more repetition
z_link_func = @(x) 1./(1+exp(x))';
trvars = dat.Properties.VariableNames;
zcells = find(~cellfun(@isempty, strfind(trvars, 'z_')));
for z = zcells,
    dat.(trvars{z}) = z_link_func(dat.(trvars{z}))';
end

% recode some stuff
params = {'v', 'z'};
for p = 1:length(params),
    tp = {'alternating', 'neutral', 'repetitive'};
    for t = 1:length(tp),
        dat.(sprintf('%s_prevcorrect_%s__regressdczprevrespstim', params{p}, tp{t})) = ...
            dat.(sprintf('%s_prevresp_%s__regressdczprevrespstim', params{p}, tp{t})) + ...
            dat.(sprintf('%s_prevstim_%s__regressdczprevrespstim', params{p}, tp{t}));
        dat.(sprintf('%s_preverror_%s__regressdczprevrespstim', params{p}, tp{t})) = ...
            dat.(sprintf('%s_prevresp_%s__regressdczprevrespstim', params{p}, tp{t})) - ...
            dat.(sprintf('%s_prevstim_%s__regressdczprevrespstim', params{p}, tp{t}));
    end
end

params = {'v_prevresp', 'v_prevstim', 'z_prevresp', 'z_prevstim', ...
    'v_prevcorrect', 'v_preverror', 'z_prevcorrect', 'z_preverror' };
for p = 1:length(params),
    
    subplot(4,4,p);
    allaccuracy     = [splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.2), findgroups(dat.subjnr(dat.transitionprob == 0.2))),  ...
        splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.5), findgroups(dat.subjnr(dat.transitionprob == 0.5))) ...
        splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.8), findgroups(dat.subjnr(dat.transitionprob == 0.8)))];
    
    allserialbias   = dat{dat.session == 0, {sprintf('%s_alternating__regressdczprevrespstim', params{p}) ...
        sprintf('%s_neutral__regressdczprevrespstim', params{p}) ...
        sprintf('%s_repetitive__regressdczprevrespstim', params{p})}};
    
    allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
    g = gscatter(allserialbias(:),  allaccuracy(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
    xlabel(regexprep(params{p}, '\_', ' ~ ')); ylabel('p(repeat)');
    axis tight;
    lsline;
    hline(0.5);
    
    % vline(0);
    
    axisNotSoTight; axis square; box off;
    
    %     % make correlations
    %     for i = 1:3,
    %         [rho(i), pval(i)] = corr(allaccuracy(:, i), allserialbias(:, i), 'type', 'spearman');
    %     end
    
end
l = legend({'Repetitive', 'Neutral', 'Alternating'});
l.Position(2) = l.Position(2) - 0.25;
l.Box = 'off';
print(gcf, '-dpdf', '~/Data/serialHDDM/fig7_Anke_conditions2.pdf');

% ======================================================= %
% ANKE'S DATA - correlate with stimcoding model
% ======================================================= %

close all; clc;
params = {'z', 'dc'};

for p = 1:length(params),
    
    subplot(4,4,p);
    allaccuracy     = [splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.2), ....
        findgroups(dat.subjnr(dat.transitionprob == 0.2))),  ...
        splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.5), ...
        findgroups(dat.subjnr(dat.transitionprob == 0.5))) ...
        splitapply(@nanmean, dat.repetition(dat.transitionprob == 0.8), ...
        findgroups(dat.subjnr(dat.transitionprob == 0.8)))];
    
    allserialbias  = dat{dat.session == 0, ...
        {sprintf('%s_1_alternating__stimcodingdczprevresp', params{p}), ...
        sprintf('%s_2_alternating__stimcodingdczprevresp', params{p}), ...
        sprintf('%s_1_neutral__stimcodingdczprevresp', params{p}), ...
        sprintf('%s_2_neutral__stimcodingdczprevresp', params{p}), ...
        sprintf('%s_1_repetitive__stimcodingdczprevresp', params{p}), ...
        sprintf('%s_2_repetitive__stimcodingdczprevresp', params{p})}};
    allserialbias = [allserialbias(:, 1) - allserialbias(:, 2) ...
        allserialbias(:, 3) - allserialbias(:, 4) ...
        allserialbias(:, 5) - allserialbias(:, 6)];
    
    allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
    g = gscatter(allserialbias(:),  allaccuracy(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
    xlabel(regexprep(params{p}, '\_', ' ~ ')); ylabel('p(repeat)');
    axis tight;
    lsline;
    hline(0.5);
    axisNotSoTight; axis square; box off;
    
end
l = legend({'Repetitive', 'Neutral', 'Alternating'});
l.Position(2) = l.Position(2) - 0.25;
l.Box = 'off';
print(gcf, '-dpdf', '~/Data/serialHDDM/fig7_Anke_conditions2_stimcoding.pdf');


% ======================================================= %
% ANKE'S DATA - correlate accuracy with serial bias
% ======================================================= %

close all; clear; clc;
dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
    'Anke_2afc_sequential'));

% get the sessions separately per condition
subplot(4,4,1);
allaccuracy     = [splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.2), findgroups(dat.subjnr(dat.transitionprob == 0.2))),  ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.5), findgroups(dat.subjnr(dat.transitionprob == 0.5))) ...
    splitapply(@nanmean, dat.accuracy(dat.transitionprob == 0.8), findgroups(dat.subjnr(dat.transitionprob == 0.8)))];

allserialbias   = dat{dat.session == 0, {'v_prevresp_alternating__regressdcprevresp', ...
    'v_prevresp_neutral__regressdcprevresp', ...
    'v_prevresp_repetitive__regressdcprevresp'}};

allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
g = gscatter(allserialbias(:), allaccuracy(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
xlabel('v ~ prevresp'); ylabel('accuracy (%)');
lsline; axisNotSoTight; axis square; box off;
xlim([-1 1]); ylim([0.6 1]); set(gca, 'ytick', 0.5:0.1:1);

% make correlations
for i = 1:3,
    [rho(i), pval(i)] = corr(allaccuracy(:, i), allserialbias(:, i), 'type', 'spearman');
end
l = legend(g, {sprintf('Repetitive, \\rho = %.2f, p = %.3f', rho(1), pval(1)), ...
    sprintf('Neutral, \\rho = %.2f, p = %.3f', rho(2), pval(2)),...
    sprintf('Alternating, \\rho = %.2f, p = %.3f', rho(3), pval(3))});
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';
title('Anke biased');

% ======================================================= %
% ANKE'S DATA - serial bias vs RT modulation
% ======================================================= %

% get the sessions separately per condition
subplot(4,4,5);
allserialbias   = dat{dat.session == 0, {'v_prevresp_alternating__regressdcprevrespstimrtpupil', ...
    'v_prevresp_neutral__regressdcprevrespstimrtpupil', ...
    'v_prevresp_repetitive__regressdcprevrespstimrtpupil'}};
allmodulation   = dat{dat.session == 0, {'v_prevrespprevrt_alternating__regressdcprevrespstimrtpupil', ...
    'v_prevrespprevrt_neutral__regressdcprevrespstimrtpupil', ...
    'v_prevrespprevrt_repetitive__regressdcprevrespstimrtpupil'}};

allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
g = gscatter(allserialbias(:), allmodulation(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
xlabel('v ~ prevresp'); ylabel('v ~ prevresp * prevrt');
lsline; axisNotSoTight; axis square; box off;
%xlim([-1 1]); ylim([0.6 1]);
%set(gca, 'ytick', 0.5:0.1:1);

% make correlations
for i = 1:3,
    [rho(i), pval(i)] = corr(allserialbias(:, i), allmodulation(:, i), 'type', 'spearman');
end
l = legend(g, {sprintf('Repetitive, \\rho = %.2f, p = %.3f', rho(1), pval(1)), ...
    sprintf('Neutral, \\rho = %.2f, p = %.3f', rho(2), pval(2)),...
    sprintf('Alternating, \\rho = %.2f, p = %.3f', rho(3), pval(3))});
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';

% ======================================================= %
% ANKE'S DATA - serial bias vs Pupil modulation
% ======================================================= %

% get the sessions separately per condition
subplot(4,4,9);
allserialbias   = dat{dat.session == 0, {'v_prevresp_alternating__regressdcprevrespstimrtpupil', ...
    'v_prevresp_neutral__regressdcprevrespstimrtpupil', ...
    'v_prevresp_repetitive__regressdcprevrespstimrtpupil'}};
allmodulation   = dat{dat.session == 0, {'v_prevrespprevpupil_alternating__regressdcprevrespstimrtpupil', ...
    'v_prevrespprevpupil_neutral__regressdcprevrespstimrtpupil', ...
    'v_prevrespprevpupil_repetitive__regressdcprevrespstimrtpupil'}};

allcondition = repmat([0.2 0.5 0.8], size(allserialbias, 1), 1);
g = gscatter(allserialbias(:), allmodulation(:), allcondition(:), flipud(linspecer(3)), '.', 8, 'off');
xlabel('v ~ prevresp'); ylabel('v ~ prevresp * prevpupil');
lsline; axisNotSoTight; axis square; box off;
%xlim([-1 1]); ylim([0.6 1]);
%set(gca, 'ytick', 0.5:0.1:1);

% make correlations
for i = 1:3,
    [rho(i), pval(i)] = corr(allserialbias(:, i), allmodulation(:, i), 'type', 'spearman');
end
l = legend(g, {sprintf('Repetitive, \\rho = %.2f, p = %.3f', rho(1), pval(1)), ...
    sprintf('Neutral, \\rho = %.2f, p = %.3f', rho(2), pval(2)),...
    sprintf('Alternating, \\rho = %.2f, p = %.3f', rho(3), pval(3))});
l.Position(1) = l.Position(1) + 0.25;
l.Box = 'off';

print(gcf, '-dpdf', '~/Data/serialHDDM/fig7_Anke_bias.pdf');

