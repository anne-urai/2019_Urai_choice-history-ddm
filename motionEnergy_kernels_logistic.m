function motionEnergy_kernels_logistic()
addpath('~/Desktop/code/gramm');

%{
    Psychophysical kernels were computed by averaging contrast fluctuations as a
function of the sample number of the test stimulus. We treated the reference grating
as if it had been shown at the same time of the test stimulus for the same duration.
The psychophysical kernel at time point t is then given by:
E(t) = hE(t)S i + hE(t)N i
E(t)S is the contrast fluctuation of the selected option at time t and E(t)N is the c
ontrast fluctuation of the non selected stimulus at time t. Absolute stimulus contrast
was transformed into contrast fluctuation by subtracting average generative contrast
values (e.g. QUEST threshold) for the respective trial type (0.5 + threshold or 0.5 ?
threshold). Average generative contrast values were computed in blocks of 100 trials,
corresponding to natural blocks in the experiment after which participants took a
short break. The expectation is across trials.
%}

path = '~/Data/psychophysicalKernels';
load(sprintf('%s/%s', path, 'motionEnergyData_AnkeMEG.mat'));

% ridge regression without crossval
[B,FitInfo] = lassoglm(data.motionenergy_normalized, (data.behavior.response > 0),'binomial', 'alpha', 0.01);

rng('default') % for reproducibility
[B,FitInfo] = lassoglm(data.motionenergy_normalized, (data.behavior.response > 0),'binomial',...
    'NumLambda',25,'CV',10);

% find best regularization
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show','Location','best') % show legend

% weights as a function of lambda
lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');

end