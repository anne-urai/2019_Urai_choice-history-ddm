function kostis_all_BICs

% ============================================ %
% BIC COMPARISON BETWEEN ALL MODELS
% ============================================ %

global mypath colors
colors(3, :) = mean(colors([1 2], :));

cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
colors = [colors; cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)]); colors(1:3, :); cols2([5 3 4], :)];

% entry point: DDM, DDM-collapsingBounds, ddm-stimOnsetAccumulation
% OU-stimOnsetAccumulation, OU-stimOnsetAccumulation-collapsingBounds

results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
baselineModel = results.ddmK_vanilla_bic;
% ddmK = normal DDM, with ramping
% ddmColl = collapsing bounds
% ddmDColl = collapsing bounds, accumulation from stimulus onset
% ouK = O-U, accumulation from stimulus offset
% ouD = O-U, accumulation from stimulus onset
% ouDColl = O-U with collapsing bounds, accumulation from stimulus onset

model(1).data = [results.ddmK_z_bic results.ddmK_dc_bic  results.ddmK_dcz_bic];
model(1).vanilla = results.ddmK_vanilla_bic;
model(1).name = {'1. Standard DDM' ''};
model(1).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(1).colors = [colors; mean(colors([1 2], :))];
model(1).basevalue = 0;
model(1).ylim = [-6500 500];
model(1).subplot = 1;

model(end+1).data = [ results.ddmK_dc_bic results.ddmK_rp_bic results.ddmK_rp2_bic];
model(end).vanilla = results.ddmK_vanilla_bic;
model(end).name = {'2. Standard DDM', 'dynamic v_{bias}'};
model(end).ticklabels = {'offset', 'ramp', 'offset+ramp'};
cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
model(end).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
model(end).basevalue = 0;
model(end).ylim = [-6500 500];
model(end).subplot = 2;

% model(3).data = [results.ddmColl_z_bic results.ddmColl_dc_bic results.ddmColl_dcz_bic];
% model(3).vanilla = results.ddmColl_vanilla_bic;
% model(3).name = {'c. Standard DDM', 'collapsing bounds'};
% model(3).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
% model(3).colors = [colors; mean(colors([1 2], :))];
% model(3).basevalue = -6000;
% model(3).ylim = [-6500 500];
% model(3).subplot = 3;
% 
% model(7).data = [results.ouK_sp_bic results.ouK_input_bic results.ouK_lambda_bic];
% model(7).vanilla = results.ouK_vanilla_bic;
% model(7).name = {'d. Offset O-U'};
% model(7).ticklabels = {'offset bias', 'input bias', '\lambda bias'};
% model(7).colors = [cols2([5 3 4], :)];
% model(7).basevalue = 0;
% model(7).ylim = [-6500 500];
% model(7).subplot = 4;

model(end+1).data = [results.ddmDColl_z_bic results.ddmDColl_dc_bic results.ddmDColl_dcz_bic];
model(end).vanilla = NaN;
model(end).name = {'3. Dynamic DDM', 'collapsing bounds'};
model(end).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(end).colors = [colors; mean(colors([1 2], :))];
model(end).basevalue = 0;
model(end).ylim = [-6500 500];
model(end).subplot = 3;

model(end+1).data = [results.ouD_sp_bic results.ouD_input_bic results.ouD_lambda_bic];
model(end).vanilla = results.ouD_vanilla_bic;
model(end).name = {'4. Dynamic O-U' ''};
model(end).ticklabels = {'offset bias', 'input bias', '\lambda bias'};
model(end).colors = [cols2([5 3 4], :)];
model(end).basevalue = 0;
model(end).ylim = [-6500 500];
model(end).subplot = 4;

model(end+1).data = [results.ouDColl_sp_bic results.ouDColl_input_bic results.ouDColl_lambda_bic];
model(end).vanilla = NaN;
model(end).name = {'5. Dynamic O-U' 'collapsing bounds'};
model(end).ticklabels = {'offset bias', 'input bias', '\lambda bias'};
model(end).colors = [cols2([5 3 4], :)];
model(end).basevalue = 0;
model(end).ylim = [-6500 500];
model(end).subplot = 5;

% model(8).data = [ results.ddmD_dc_bic results.ddmD_rp_bic results.ddmD_rp2_bic];
% model(8).vanilla = results.ddmD_vanilla_bic;
% model(8).name = {'g. Dynamic DDM', 'dynamic v_{bias}'};
% model(8).ticklabels = {'offset', 'ramp', 'offset+ramp'};
% cols1 = cbrewer('qual', 'Set1', 8);
% cols2 = cbrewer('qual', 'Dark2', 8);
% model(8).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
% model(8).basevalue = 0;
% model(8).ylim = [0 2000];
% model(8).subplot = 7;


% move subplots closer together
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.1 0.01], [0.1 0.01]);

close all;
for m = 1:length(model),
    
    %everything relative to the full model
    mdldic = bsxfun(@minus, mean(model(m).data), mean(baselineModel));
    subplot(2,5,model(m).subplot);     
    hold on;
    plot([0.5 length(mdldic)+0.5], ...
        [mean(model(m).vanilla)- mean(baselineModel) mean(model(m).vanilla)-mean(baselineModel)], '-.k');
    
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
    %ylim(model(m).ylim);
    ylim([-4700 500]);
    disp(get(gca, 'ylim'))
    
    if m > 2
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
                'horizontalalignment', 'center', 'color', 'k');
        end
    end
    
    box off;
    set(gca, 'color', 'none');
    set(gca, 'xcolor', 'k', 'ycolor', 'k');
    offsetAxes;
    set(gca, 'xtick', 1:length(model(m).ticklabels), 'xticklabel', model(m).ticklabels, 'xticklabelrotation', -30);
    title(model(m).name);
    axis square;
    
    switch m
        case 1
            ylabel({'\DeltaBIC from DDM'; 'without history'}, 'interpreter', 'tex');
        %case 3
          %  set(gca, 'ytick', [-6400:200:-6000], 'yticklabel', {'-6400', '-6200', '0'});
        otherwise
            set(gca, 'ycolor', 'w', 'yticklabel', []);
    end
end

%xlabel('DDM');
drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/BICs_kostis_all.pdf'));
end

