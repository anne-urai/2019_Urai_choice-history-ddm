function kostis_summary
% check Kostis' new set of results, 30 August 2018

close all; clc;
global mypath datasets datasetnames
kostisPath = '~/Data/HDDM/Anke_MEG_transition/newDDM/';
cd(kostisPath);

clear;
load allmodels_2.mat;
load history;

LL=outgf; % likelihoods for each subject
% rows are: 1-naive, 2-DC, 3-startint point
% 4-ramping with slope, 5-ramping with 2 params
% 6- both DC and starting point;

%% GRAB ALL PARAMETERS

load history;
kk(6,:)=[]; % remove this subject

summary.repetition = kk(:, 2); % not sure how this is computed
summary.vshift_dc_prevresp = params2(:,end);
summary.zshift_z_prevresp = params3(:,end);
summary.vramp_dc_prevresp = params4(:,end);
summary.vramp2_offset = params5(:,end);
summary.vramp2_slope = params5(:,end-1).*sign(params5(:,end));

summary.vshift_dc_z_prevresp = params6(:,end-1);
summary.zshift_dc_z_prevresp = params6(:,end);

close all;
corrplot(summary, {'zshift_z_prevresp', 'vshift_dc_prevresp'},  {'repetition'});
suptitle('Independent fits, zshift vs. vshift');
[rhodiff, pval] = steigers(summary.repetition, summary.zshift_z_prevresp, summary.vshift_dc_prevresp);
suplabel(sprintf('Rho difference = %.2f, p = %.3f', rhodiff, pval), 'x');
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/kostis_check_1.pdf'));

close all;
corrplot(summary, {'zshift_dc_z_prevresp', 'vshift_dc_z_prevresp'},  {'repetition'});
suptitle('Joint fits, zshift vs. vshift');
[rhodiff, pval] = steigers(summary.repetition, summary.zshift_dc_z_prevresp, summary.vshift_dc_z_prevresp);
suplabel(sprintf('Rho difference = %.2f, p = %.3f', rhodiff, pval), 'x');
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/kostis_check_2.pdf'));

close all;
corrplot(summary, {'vshift_dc_prevresp', 'vramp_dc_prevresp'},  {'repetition'});
suptitle('Independent fits, offset vs. ramping vbias');
[rhodiff, pval] = steigers(summary.repetition, summary.vshift_dc_prevresp, summary.vramp_dc_prevresp);
suplabel(sprintf('Rho difference = %.2f, p = %.3f', rhodiff, pval), 'x');
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/kostis_check_3.pdf'));

close all;
corrplot(summary, {'vramp2_offset', 'vramp2_slope'},  {'repetition'});
suptitle('Joint fits, offset vs. ramping vbias');
[rhodiff, pval] = steigers(summary.repetition, summary.vramp2_offset, summary.vramp2_slope);
suplabel(sprintf('Rho difference = %.2f, p = %.3f', rhodiff, pval), 'x');
tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/kostis_check_4.pdf'));

end

function [rhodiff, pval] = steigers(var1,var2,var3)

rho1 = corr(var2, var1, 'rows', 'complete', 'type', 'spearman');
rho2 = corr(var3, var1, 'rows', 'complete', 'type', 'spearman');
rho3 = corr(var1, var2, 'rows', 'complete', 'type', 'spearman');

[rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3, numel(var1), 0.05);

end
