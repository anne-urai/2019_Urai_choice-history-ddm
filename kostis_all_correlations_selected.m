function kostis_all_correlations_selected

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

% ============================================ %

global mypath colors
colors(3, :) = mean(colors([1 2], :));

cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
colors = [colors; cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)]); colors(1:3, :); cols2([5 3 4], :)];
grey = [0.7 0.7 0.7];

% entry point: DDM, DDM-collapsingBounds, ddm-stimOnsetAccumulation
% OU-stimOnsetAccumulation, OU-stimOnsetAccumulation-collapsingBounds

% results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));
baselineModel = results.repetitionK;

% ddmK = normal DDM, with ramping
% ddmColl = collapsing bounds
% ddmDColl = collapsing bounds, accumulation from stimulus onset
% ouK = O-U, accumulation from stimulus offset
% ouD = O-U, accumulation from stimulus onset
% ouDColl = O-U with collapsing bounds, accumulation from stimulus onset

model(1).data       = {{results.ddmK_dcz_zbias results.ddmK_dcz_dcbias}};
model(1).name       = {'1. Standard' 'DDM'};
model(1).ticklabels = {{'z', 'v_{bias}'}};
model(1).colors     = [grey; grey];
model(1).subplot    = 1;

model(end+1).data = {{results.ddmK_rp2_offset results.ddmK_rp2_slope}};
model(end).name = {'2. DDM', 'dynamic v_{bias}'};
model(end).ticklabels = {{'constant', 'ramp'}};
cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
model(end).colors = [grey; grey];
model(end).subplot = 2;

% model(end+1).data = {{results.ddmD_rp2_offset results.ddmD_rp2_slope}};
% model(end).name = {'2b. Dynamic DDM', 'drift bias'};
% model(end).ticklabels = {{'constant', 'ramp'}};
% cols1 = cbrewer('qual', 'Set1', 8);
% cols2 = cbrewer('qual', 'Dark2', 8);
% model(end).colors = [grey; grey];
% model(end).subplot = 3;

% model(end+1).data = {results.ddmColl_z_zbias results.ddmColl_dc_dcbias {results.ddmColl_dcz_zbias results.ddmColl_dc_dcbias}};
% model(3).vanilla = results.ddmColl_vanilla_bic;
% model(3).name = {'c. Standard DDM', 'collapsing bounds'};
% model(3).ticklabels = {'z', 'v_{bias}', {'z_{bias}', 'v_{bias}'}};
% model(3).colors = [colors; mean(colors([1 2], :))];
% model(3).subplot = 3;
%
% model(7).data = {results.ouK_sp_spbias results.ouK_input_inputbias results.ouK_lambda_lambdabias};
% model(7).vanilla = results.ouK_vanilla_bic;
% model(7).name = {'d. Offset O-U'};
% model(7).ticklabels = {'offset bias', 'input bias', '\lambda bias'};
% model(7).colors = [cols2([5 3 4], :)];
% model(7).subplot = 4;

model(end+1).data = {{results.ddmDColl_dcz_zbias results.ddmDColl_dc_dcbias}};
model(end).name = {'3. Dynamic DDM', 'collapsing bounds'};
model(end).ticklabels = {{'z', 'v_{bias}'}};
model(end).colors = [grey; grey];
model(end).subplot = 3;

% model(end+1).data = {results.ouD_lambda_lambdabias};
% model(end).vanilla = NaN;
% model(end).name = {'4. Dynamic O-U' ''};
% model(end).ticklabels = {'\lambda bias'};
% model(end).colors = grey;
% model(end).subplot = 4;

model(end+1).data = {results.ouDColl_input_inputbias};
model(end).vanilla = NaN;
model(end).name = {'4. Leaky accumulator' 'collapsing bounds'};
model(end).ticklabels = {'input bias'};
model(end).colors = grey;
model(end).subplot = 4;

% model(8).data = {results.ddmD_dc_dcbias	results.ddmD_rp_slope {results.ddmD_rp2_offset results.ddmD_rp2_slope}};
% model(8).name = {'g. Dynamic DDM', 'dynamic v_{bias}'};
% model(8).ticklabels = {'offset', 'ramp', {'offset', 'ramp'}};
% cols1 = cbrewer('qual', 'Set1', 8);
% cols2 = cbrewer('qual', 'Dark2', 8);
% model(8).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
% model(8).subplot = 7;

% move subplots closer together
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.1 0.01], [0.1 0.01]);

close all;
set(gcf, 'defaultaxesfontsize', 4, 'defaultaxestitlefontsizemultiplier', 1);

for m = 1:length(model),
    
    %everything relative to the full model
    subplot(4,max([model(:).subplot])+3,model(m).subplot);
    hold on;
    xticks = [];
    xticklabels = {};
    
    %%%%%%
    for i = 1:length(model(m).data),
        if ~iscell(model(m).data{i})
            [rho, ci, pval, bf] = spearmans(model(m).data{i}, baselineModel);
            b = bar(i, rho, 'facecolor', model(m).colors(i, :), 'barwidth', 0.6, 'BaseValue', 0, ...
                'edgecolor', 'none');
            errorbar(i, rho, ci(1)-rho, ci(2)-rho, 'k', 'marker', 'none','capsize', 0, 'linewidth', 1);
            mysigstar(gca, i, 0.05, pval, 'w');
            xticks = [xticks i];
            xticklabels = [xticklabels model(m).ticklabels{i}];
            
        else
            % for the two separately
            xshift = 0.25;
            [rho1, ci, pval, bf] = spearmans(model(m).data{i}{1}, baselineModel);
            b = bar(i-xshift, rho1, 'facecolor', model(m).colors(i, :), 'barwidth', 0.45, 'BaseValue', 0, ...
                'edgecolor', 'none');
            errorbar(i-xshift, rho1, ci(1)-rho1, ci(2)-rho1, 'k', 'marker', 'none','capsize', 0, 'linewidth', 1);
            mysigstar(gca, i-xshift, 0.05, pval, 'w');
            xticks = [xticks i-xshift];
            xticklabels = [xticklabels model(m).ticklabels{i}{1}];
            pval1 = pval;
            
            % second one
            [rho2, ci, pval, bf] = spearmans(model(m).data{i}{2}, baselineModel);
            b = bar(i+xshift, rho2, 'facecolor', model(m).colors(i, :), 'barwidth', 0.45, 'BaseValue', 0, ...
                'edgecolor', 'none');
            errorbar(i+xshift, rho2, ci(1)-rho2, ci(2)-rho2, 'k', 'marker', 'none','capsize', 0, 'linewidth', 1);
            mysigstar(gca, i+xshift, 0.05, pval, 'w');
            xticks = [xticks i+xshift];
            xticklabels = [xticklabels model(m).ticklabels{i}{2}];
            pval2 = pval;
            
            % difference between them
            % compute the difference in correlation
            [rho3] = corr(model(m).data{i}{1}, model(m).data{i}{2}, ...
                'rows', 'complete', 'type', 'spearman');
            [~, ~, pval] = rddiffci(rho1,rho2,rho3, length(baselineModel), 0.05);
            mysigstar(gca, [i-xshift-0.1 i+xshift+0.1], -0.5, pval);
            
            disp(sprintf('%s, %s rho %f p %s, %s rho %f p %f', ...
                model(m).name{1}, model(m).ticklabels{i}{1}, rho1, pval1, model(m).ticklabels{i}{2}, rho2, pval2));
            
        end
    end
    
    % now break the y-axes
    xlim([0.5 length(model(m).data)+0.5]);
    hline(0, [0 0 0]);
    ylim([-0.5 1]);
    box off;
    set(gca, 'color', 'none');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    offsetAxes;
    set(gca, 'xtick', xticks, 'xticklabel', xticklabels, 'xticklabelrotation', -30);
    title(model(m).name, 'fontweight', 'normal', 'fontangle', 'italic', 'fontsize', 8);
    %axis square;
    
    switch model(m).subplot
        case {1, 7}
            ylabel({'Correlation \rho'; 'with P(repeat)'}, 'interpreter', 'tex');
        otherwise
            set(gca, 'ycolor', 'w', 'yticklabel', []);
    end
end

drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/Correlations_kostis_selected.pdf'));
end

