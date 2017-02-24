clear all; close all;
usepath = '~/Data/RT_RDK/HDDM';
% mdls = {'stimcoding', 'regress_dc'};

% ============================================ %
% GET PARAMETERS OF INTEREST
% ============================================ %

stimcoding  = readtable(sprintf('%s/%s/results-md1.csv', usepath, 'stimcoding'), 'readrownames', 1);
regressdc  = readtable(sprintf('%s/%s/results-md1.csv', usepath, 'regress_dc'), 'readrownames', 1);

params = {'a_subj', 't_subj', 'z_subj'}

% boundary separation
for p = 1:length(params),
    subplot(3,3,p);
    plot(stimcoding{strncmp(stimcoding.Properties.RowNames, params{p}, numel(params{p})), 2}, ...
        regressdc{strncmp(regressdc.Properties.RowNames,  params{p}, numel(params{p})), 2}, '.');
    axisNotSoTight; lsline; box off; %axis square;
    title(params{p}, 'interpreter', 'none');
end

suplabel('stimcoding', 'x'); suplabel('regression', 'y');

 saveas(gcf, sprintf('%s/Figures/HDDMcheck.pdf', usepath), 'pdf');
