function kostis_all_BICs

% ============================================ %
% BIC COMPARISON BETWEEN ALL MODELS
% ============================================ %

global mypath colors
colors(3, :) = mean(colors([1 2], :));

cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
colors = [colors; cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)]); colors(1:3, :); cols2([5 3 4], :)];

results = readtable('/Users/urai/Data/HDDM/summary/Anke_MEG_transition/allindividualresults_kostis.csv');
baselineModel = results.ddmK_vanilla_bic;

model(1).data = [results.ddmK_z_bic results.ddmK_dc_bic  results.ddmK_dcz_bic];
model(1).name = {'DDM', 'standard'};
model(1).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(1).colors = [colors; mean(colors([1 2], :))];
model(1).basevalue = 0;
model(1).ylim = [-3000 300];

model(2).data = [ results.ddmK_dc_bic results.ddmK_rp_bic results.ddmK_rp2_bic];
model(2).name = {'DDM', 'dynamic v_{bias}'};
model(2).ticklabels = {'offset', 'ramp', 'offset+ramp'};
cols1 = cbrewer('qual', 'Set1', 8);
cols2 = cbrewer('qual', 'Dark2', 8);
model(2).colors = [cols1(2, :); cols2(6, :); nanmean([cols1(2, :); cols2(6, :)])];
model(2).basevalue = 0;
model(2).ylim = [-3000 300];

model(3).data = [results.ddmColl_z_bic results.ddmColl_dc_bic results.ddmColl_dcz_bic];
model(3).name = {'DDM', 'collapsing bounds'};
model(3).ticklabels = {'z', 'v_{bias}', 'z+v_{bias}'};
model(3).colors = [colors; mean(colors([1 2], :))];
model(3).basevalue = -6000;
model(3).ylim = [-6400 -5950];

model(4).data = [results.ouD_input_bic results.ouD_lambda_bic];
model(4).name = {'Ornstein-' 'Uhlenbeck'};
model(4).ticklabels = {'input bias', '\lambda bias'};
model(4).colors = [cols2([3 4], :)];
model(4).basevalue = 0;
model(4).ylim = [-3000 300];

% move subplots closer together
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.1 0.01], [0.1 0.01]);

close all;
for m = 1:length(model),
    
    %everything relative to the full model
    mdldic = bsxfun(@minus, mean(model(m).data), mean(baselineModel));
    subplot(4,6,m);
    
    %%%%%%
    hold on;
    for i = 1:length(mdldic),
        b = bar(i, mdldic(i), 'facecolor', model(m).colors(i, :), 'barwidth', 0.7, 'BaseValue', model(m).basevalue, ...
            'edgecolor', 'none');
    end

    % now break the y-axes
    xlim([0.5 length(mdldic)+0.5]);
    hline(model(m).basevalue, [0 0 0]);
    ylim(model(m).ylim);
    
    if ismember(m, [3 4]),
        for i = 1:length(mdldic)
            text(i, min([0 mdldic(i)]) + 0.1*range(get(gca, 'ylim')), ...
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
    
    switch m
        case 1
            ylabel({'\DeltaBIC from DDM'; 'without history'}, 'interpreter', 'tex');
        case 3
            set(gca, 'ytick', [-6400:200:-6000], 'yticklabel', {'-6400', '-6200', '0'});
        otherwise
            set(gca, 'ycolor', 'w', 'yticklabel', []);
    end
end

%xlabel('DDM');
drawnow; tightfig;
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/BICs_kostis_all.pdf'));
end

