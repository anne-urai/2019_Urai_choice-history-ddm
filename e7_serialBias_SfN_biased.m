% ======================================================= %
% ANKE'S DATA - correlate accuracy with serial bias
% ======================================================= %
close all; clear;
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

