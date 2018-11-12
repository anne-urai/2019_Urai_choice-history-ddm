function meg_regression_posteriors

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global mypath

datasets        = {'MEG_MEGdata'};
d = 1;
mdls = {'regress_dc_z_motorstart', ...
    'regress_dc_z_visualgamma'};

for m = 1:length(mdls),
    close all;
    dat = readtable(sprintf('%s/%s/%s/group_traces.csv', mypath, datasets{d}, mdls{m}));
    
    if strfind(mdls{m}, 'motorstart'),
        dat.v_motorstart = dat.v_motorbeta;
        dat.z_motorstart = dat.z_motorbeta;
    end
    
    
    % which parameters do we need?
    params = regexprep(mdls{m}, 'regress_dc_z_', '');
    params = strsplit(params, '_');
    for p = 1:length(params),
        
        switch params{p}
            case 'visualgamma'
                paramname = 'visual \gamma';
            case 'motorstart'
                paramname = 'motor \beta baseline';
        end
        
        sp1 = subplot(4,4,(p-1)*4+1); hold on;
        plotHist(dat.(['v_' params{p}]));
        xlabel(['v_{bias} ~ ' paramname]);
        
        set(gca, 'yticklabel', []);
        ylabel({'Posterior' 'probability'});
        % second one
        subplot(4,4,(p-1)*4+2); hold on;
        plotHist(dat.(['z_' params{p}]));
        xlabel(['z ~ ' paramname]);
        set(gca, 'yticklabel', []);
        
    end
    
    sp1.Position(1) = sp1.Position(1) + 0.03;
    % suplabel('Posterior probability', 'y');
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/MEG_posteriors_%s.pdf', mdls{m}));
end

end

function plotHist(y, param)

y = y*10^3;
% make a nice-looking plot
[xHist,yHist] = ksdensity(y);

% find the percentiles and make patch
xPatch = xHist;
yPatch = yHist;

lowerbnd = prctile(y, 2.5);
upperbnd = prctile(y, 97.5);
yPatch(yPatch < lowerbnd) = lowerbnd;
yPatch(yPatch > upperbnd) = upperbnd;

patch(yPatch, xPatch, [0.7 0.7 0.7], ...
    'facecolor', [0.7 0.7 0.7], 'edgecolor', 'none', 'facealpha', 0.5);
hold on;

% outline mean
plot([mean(y) mean(y)], [min(xHist) max(xHist)], 'w-', 'linewidth', 2);

% outline thicker, on top
plot(yHist, xHist, '-', 'color', [0.3 0.3 0.3], 'linewidth', 2);

% hold on; histogram(dat.(param), 'edgecolor', 'none', 'facecolor', ...
%     [0.5 0.5 0.5]);

red = linspecer(2);
axis tight;
vline(0, 'color', red(2, :));
offsetAxes;
text(mean([mean(get(gca, 'xlim')) 0.75*max(get(gca, 'xlim'))]), mean(get(gca, 'ylim')), ...
    sprintf('p = %.3f', min([mean(y < 0) mean(y > 0)])), 'fontsize', 6);


end