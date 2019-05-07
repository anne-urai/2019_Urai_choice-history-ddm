function kostis_all_AICs

% ============================================ %
% AIC COMPARISON BETWEEN ALL MODELS
% ============================================ %

global mypath colors

cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
colors = [colors; cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)]); colors(1:3, :); cols2([5 3 4], :)];

% entry point: DDM, DDM-collapsingBounds, ddm-stimOnsetAccumulation
% OU-stimOnsetAccumulation, OU-stimOnsetAccumulation-collapsingBounds

% results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
results = readtable(sprintf('%s/summary/%s/allindividualresults_kostis.csv', mypath, 'Anke_MEG_transition'));
% baselineModel = results.ddmK_vanilla_aic;

% ddmK = normal DDM, with ramping
% ddmColl = collapsing bounds
% ddmDColl = collapsing bounds, accumulation from stimulus onset
% ouK = O-U, accumulation from stimulus offset
% ouD = O-U, accumulation from stimulus onset
% ouDColl = O-U with collapsing bounds, accumulation from stimulus onset

model(1).data = [results.ddmK_z_aic results.ddmK_dc_aic  results.ddmK_dcz_aic];
model(1).vanilla = results.ddmK_vanilla_aic;
model(1).name = {'1. Standard' 'DDM'};
model(1).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(1).colors = [colors; mean(colors([1 2], :))];
model(1).basevalue = unique(results.ddmK_vanilla_aic);
model(1).ylim = [-650 0];
model(1).subplot = 1;

model(end+1).data = [ results.ddmK_dc_aic results.ddmK_rp_aic results.ddmK_rp2_aic];
model(end).vanilla = results.ddmK_vanilla_aic;
model(end).name = {'2. DDM', 'dynamic v_{bias}'};
model(end).ticklabels = {'constant v_{bias}', 'ramp', 'constant+ramp'};
cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
model(end).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
model(end).basevalue = unique(results.ddmK_vanilla_aic);
model(end).ylim = [-650 0];
model(end).subplot = 2;

% model(end+1).data = [ results.ddmD_dc_aic results.ddmD_rp_aic results.ddmD_rp2_aic];
% model(end).vanilla = results.ddmD_vanilla_aic;
% model(end).name = {'2b. Dynamic DDM', 'drift bias'};
% model(end).ticklabels = {'constant v_{bias}', 'ramp', 'constant+ramp'};
% cols1 = cbrewer('qual', 'Set1', 8);
% cols2 = cbrewer('qual', 'Dark2', 8);
% model(end).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
% model(end).basevalue = unique(results.ddmD_vanilla_aic);
% model(end).ylim = [-8000 0];
% model(end).subplot = 4;

model(end+1).data = [results.ddmDColl_z_aic results.ddmDColl_dc_aic results.ddmDColl_dcz_aic];
model(end).vanilla = unique(results.ddmD_vanilla_aic)-6000;
model(end).name = {'3. Dynamic DDM', 'collapsing bounds'};
model(end).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(end).colors = [colors; mean(colors([1 2], :))];
model(end).basevalue = unique(results.ddmD_vanilla_aic);
model(end).ylim = [-8000 -6000];
model(end).subplot = 4;

% model(end+1).data = [results.ouD_sp_aic results.ouD_input_aic results.ouD_lambda_aic];
% model(end).vanilla = results.ddmK_vanilla_aic;
% model(end).name = {'4. Leaky' 'accumulator'};
% model(end).ticklabels = {'starting point bias', 'input bias', '\lambda bias'};
% model(end).colors = [cols2([5 3 4], :)];
% model(end).basevalue = unique(results.ddmD_vanilla_aic);
% model(end).ylim = [-8000 0];
% model(end).subplot = 5;

model(end+1).data = [results.ouDColl_sp_aic results.ouDColl_input_aic results.ouDColl_lambda_aic];
model(end).vanilla = unique(results.ddmD_vanilla_aic)-6000;
model(end).name = {'4. Leaky accumulator' 'collapsing bounds'};
model(end).ticklabels = {'starting point bias', 'input bias', '\lambda bias'};
model(end).colors = [cols2([5 3 4], :)];
model(end).basevalue = unique(results.ddmD_vanilla_aic);
model(end).ylim = [-8000 -6000];
model(end).subplot = 5;

% move subplots closer together
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.1 0.01], [0.1 0.01]);

close all;
set(gcf, 'defaultaxesfontsize', 4, 'defaultaxestitlefontsizemultiplier', 1);

for m = 1:length(model),
    
    %everything relative to the full model
    mdldic = bsxfun(@minus, mean(model(m).data), mean(model(m).basevalue));
    subplot(2,max([model(:).subplot]), model(m).subplot);     
    hold on;
    
    vn = unique(model(m).vanilla) - model(m).basevalue;
    plot([0.5 3.5], [vn vn], 'k-');
    
    %%%%%%
    for i = 1:length(mdldic),
        b = bar(i, mdldic(i), 'facecolor', model(m).colors(i, :), 'barwidth', 0.7, 'BaseValue', 0, ...
            'edgecolor', 'none');
    end
    % outline the best fit
    [~, i] = min(mdldic);
    bar(i, mdldic(i), 'facecolor', model(m).colors(i, :), 'barwidth', 0.7, 'BaseValue', 0, ...
        'edgecolor', 'k');
    
    % now break the y-axes
    xlim([0.5 length(mdldic)+0.5]);
    hline(0, [0 0 0]);
    ylim(model(m).ylim);
    disp(get(gca, 'ylim'))
    
    if m > 0
        for i = 1:length(mdldic)
            text(i, min([0 mdldic(i)]) + 0.05*range(get(gca, 'ylim')), ...
                num2str(round(mdldic(i))), ...
                'VerticalAlignment', 'bottom', 'FontSize', 4, ...
                'horizontalalignment', 'center', 'color', 'w');
        end
    else
        for i = 1:length(mdldic)
            text(i, min([0 mdldic(i)]) - 0.1*range(get(gca, 'ylim')), ...
                num2str(round(mdldic(i))), ...
                'VerticalAlignment', 'bottom', 'FontSize', 4, ...
                'horizontalalignment', 'center', 'color', 'w');
        end
    end
    
    box off;
    set(gca, 'color', 'none');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    offsetAxes;
    set(gca, 'xtick', 1:length(model(m).ticklabels), ...
        'xticklabel', model(m).ticklabels, 'xticklabelrotation', -30);
    title(model(m).name, 'fontweight', 'normal', 'fontangle', 'italic');
    axis square;
    
    switch m
        case 1
            ylabel({'\DeltaAIC from DDM'; 'without history'}, 'interpreter', 'tex');
        case 3
            ylabel({'\DeltaAIC from dyn-DDM'; 'without history'}, 'interpreter', 'tex');
        otherwise
            set(gca, 'ycolor', 'w', 'yticklabel', []);
    end
end

%xlabel('DDM');
drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/AIC_kostis_all.pdf'));
end

