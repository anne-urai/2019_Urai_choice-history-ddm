function e6_serialBias_SfN_modelFree_Reciprobit

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
datasetnames = {'RT', '2IFC', 'NatComm', 'Anke neutral'};

set(groot, 'defaultaxesfontsize', 8, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'normal', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ========================================== %
% MODELFREE MEASURE OF BIAS
% RT DISTRIBUTIONS
% ========================================== %

for d = 1:length(datasets),
    clearvars -except d datasets cnt datasetnames whichSJ ylims;
    
    % load individual serial choice bias
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    serialbias = dat.v_prevresp__regressdcprevrespstim(dat.session == 0);
    
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    csvfile = csvfile(arrayfun(@(x) ~strcmp(x.name(1),'.'), csvfile)); % remove hidden files
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    
    if d > 2, % in case of natcomm and anke's datasets
        alldata.stimulus = sign(alldata.motionenergy);
    end
    
    % add repeat
    alldata.repeat = [NaN; (diff(alldata.response) == 0)];
    stim = alldata.stimulus; alldata.stimulus(alldata.stimulus == -1) = 0;
    alldata.correct = (alldata.response == stim);
    
    subjects = unique(alldata.subj_idx)';
    [~, idx] = sort(serialbias); idx = idx';
    for sj = idx,
        
        clf; subplot(2,2,1);
        
        % from http://www.neural-code.com/index.php/tutorials/action/reaction-time/83-reciprobit-distribution
        data = alldata(alldata.subj_idx == subjects(sj) ...
            & alldata.correct == 1 & alldata.prevstim == alldata.prevresp, :);
        plotRTdistributions(data);
        
        % layout and save
        print(gcf, '-depsc', sprintf('~/Data/serialHDDM/fig4_reciprobit_%s_P%02d.eps', ...
            datasetnames{d}, sj));
    end
end
end

function h = plotRTdistributions(data)

h = [];
colors = cbrewer('div', 'PuOr', 5);
colors = colors([1 end], :);
hold on;
data.rt = data.rt * 1000;
prevresps = [0 1];

% do separate fits for choice == 1 and choice == -1 (error and correct,
% depending on the stimulus)
cnt = 1;
for p = 1:length(prevresps),
    thisdat = data(data.repeat == prevresps(p), :);
    h{cnt} = reciprobit(thisdat.rt, colors(cnt, :));
    cnt = cnt + 1;
end

end

function [h, b] = reciprobit(rt, color)

rtinv = 1./ rt;
% raw data
x = -1./sort((rt)); % multiply by -1 to mirror abscissa
n = numel(rtinv); % number of data points
y = pa_probit((1:n)./n); % cumulative probability for every data point converted to probit scale
% plot(x,y,'k.');
hold on

% quantiles
p		= [1 2 5 10:10:90 95 98 99]/100;
probit	= pa_probit(p);
q		= quantile(rt,p);
q		= -1./q;
xtick	= sort(-1./(500+[0 pa_oct2bw(100,-1:5)])); % some arbitrary xticks

% then datapoints on top
plot(q,probit,'o','markeredgecolor', color, 'markerfacecolor', 'w', 'markersize', 3);
set(gca,'XTick',xtick,'XTickLabel',-1./xtick);
xlim([min(xtick) max(xtick)]);
set(gca,'YTick',probit,'YTickLabel',p*100);
ylim([pa_probit(0.1/100) pa_probit(99.9/100)]);
axis square;
box off
xlabel('Reaction time (ms)');
ylabel('Cumulative probability');
title('Probit ordinate');

% this should be a straight line
x = q;
y = probit;
b = regstats(y,x);
h = pa_regline(b.beta,'k-');
set(h,'Color', color,'LineWidth',1);

end

function chi = pa_probit(cdf)

% CHI = PA_PROBIT(CDF)
%
% The probit function is the quantile function, i.e., the inverse
% cumulative distribution function (CDF), associated with the standard
% normal distribution.
%
% This is useful for plotting reaction times.
%

% 2013 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

myerf       = 2*cdf - 1;
myerfinv    = sqrt(2)*erfinv(myerf);
chi         = myerfinv;
end


function h = pa_regline(beta, style)
% PA_REGLINE(BETA)
%
% Plot regression line with parameters BETA through current axis.
%
% PA_REGLINE(...,'LineSpec') uses the color and linestyle specified by
% the string 'LineSpec'. See PLOT for possibilities.
%
% H = PA_REGLINE(...) returns a vector of lineseries handles in H.
%
% For example:
%
% X		= rand(100,1);
% Y		= 2*X+3+randn(size(X));
% b		= regstats(Y,X);
% beta	= b.beta;
%
% plot(X,Y,'k.');
% pa_regline(beta);
%
%
% See also REFLINE, PA_HORLINE, PA_VERLINE, PA_UNITYLINE
%

% (c) 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com


%% Initialization
if nargin < 2,
    style = 'k--';
end
if nargin < 1,
    beta = [0 1]';
end

%% Axis
x_lim	= get(gca,'Xlim');
oldhold = ishold;
hold on

%% Data
X = [ones(size(x_lim))' x_lim'];
Y = X*beta;
h       = plot(x_lim, Y, style);

%% Return
if ~oldhold
    hold off;
end

end

function F = pa_oct2bw(F1,oct)
% F2 = PA_OCT2BW(F1,OCT)
%
% Determine frequency F2 that lies OCT octaves from frequency F1
%

% (c) 2011-05-06 Marc van Wanrooij
F = F1 .* 2.^oct;
end
