
usepath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/HDDM';
close all;
fz = 7;
set(groot, 'defaultaxesfontsize', fz, ...
    'defaultaxestitlefontsizemultiplier', 1, 'defaultaxestitlefontweight', 'normal');

mdlz  = readtable(sprintf('%s/%s/group_traces.csv', usepath, 'stimcoding_prevresp_z'));

% transform back to 0.5
back2log = @(x) exp(x) ./ (1 + exp(x));
mdlz.z_trans_1_0_ = back2log(mdlz.z_trans_1_0_);
mdlz.z_trans__1_0_ = back2log(mdlz.z_trans__1_0_);

mdldc  = readtable(sprintf('%s/%s/group_traces.csv', usepath, 'stimcoding_prevresp_dc'));

% first, let z change with previous response - sample two separate posteriors
subplot(441); hold on;
histogram(mdlz.z_trans_1_0_, 'displaystyle', 'stairs');
histogram(mdlz.z_trans__1_0_, 'displaystyle', 'stairs');
xlabel('Starting point (z)'); 
title(sprintf('p = %.3f', mean(mdlz.z_trans_1_0_ < mdlz.z_trans__1_0_)));
axis tight; axis square; set(gca, 'yticklabel', []); 
ylabel('Separate posteriors');

subplot(442); hold on;
histogram(mdldc.dc_1_0_, 'displaystyle', 'stairs');
histogram(mdldc.dc__1_0_, 'displaystyle', 'stairs');
xlabel('Drift criterion (dc)');
title(sprintf('p = %.3f', mean(mdldc.dc_1_0_ < mdldc.dc__1_0_)));
axis tight; axis square; set(gca, 'yticklabel', []); 
l = legend('Previous choice 1', 'Previous choice -1'); legend boxoff;
l.Position(1) = l.Position(1) + 0.2;

% % then, look at regression models
% regress_dc  = readtable(sprintf('%s/%s/group_traces.csv', usepath, 'regress_dc_prevresp'));
% 
% subplot(444);
% histogram(regress_dc.v_prevresp, 'displaystyle', 'stairs');
% title(sprintf('p = %0.3f', mean(regress_dc.v_prevresp < 0)));
% xlabel('Previous choice -> dc'); box off; 
% axis tight; axis square; set(gca, 'yticklabel', []); 
% ylabel('Regression models');

% then, look at regression models
regress_dc  = readtable(sprintf('%s/%s/group_traces.csv', usepath, 'regress_dc_prevresp_prevpupil_prevrt'));

subplot(449);
histogram(regress_dc.v_prevresp, 'displaystyle', 'stairs');
title(sprintf('p = %0.3f', mean(regress_dc.v_prevresp < 0)));
xlabel('prevChoice'); box off; 
axis tight; axis square; set(gca, 'yticklabel', []); 
ylabel('Regression on DC');

subplot(4,4,10);
histogram(regress_dc.v_prevresp_prevpupil, 'displaystyle', 'stairs');
title(sprintf('p = %0.3f', mean(regress_dc.v_prevresp_prevpupil > 0)));
xlabel('prevChoice x prevPupil'); box off; 
axis tight; axis square; set(gca, 'yticklabel', []); 

subplot(4,4,11);
histogram(regress_dc.v_prevresp_prevrt, 'displaystyle', 'stairs');
title(sprintf('p = %0.3f', mean(regress_dc.v_prevresp_prevrt > 0)));
xlabel('prevChoice x prevRT'); box off;
axis tight; axis square; set(gca, 'yticklabel', []); 

subjectdata = subjectspecifics('GA');
print(gcf, '-dpdf', sprintf('%s/HDDM_results.pdf', subjectdata.figsdir));

% see if there is a correlation between the RT or pupil interaction weights
% and the overall weight

% clearvars -except usepath
% subplot(449);
% regress_dc      = readtable(sprintf('%s/%s/results-combined.csv', usepath, 'regress_dc_prevresp_prevpupil_prevrt'));
% prevchoiceIdx   = find(~cellfun(@isempty, strfind(regress_dc{:, 1}, 'v_prevresp_subj')));
% prevRTIdx       = find(~cellfun(@isempty, strfind(regress_dc{:, 1}, 'v_prevresp:prevrt_subj')));
% prevPupilIdx    = find(~cellfun(@isempty, strfind(regress_dc{:, 1}, 'v_prevresp:prevpupil_subj')));

