function e1b_serialBias_SfN_ModelFreeCorrelation_grey
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = length(datasets):-1:1
    disp(datasets{d});
    
    colors = linspecer(5); % red blue green
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    allresults = struct(); alltitles = {};
    switch datasets{d}
        case {'Anke_2afc_sequential'}
            
            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY
            
            % use the stimcoding difference only from alternating
            allresults(1).z_prevresp        = results.z_1_neu__stimcodingdczprevresp - results.z_2_neu__stimcodingdczprevresp;
            allresults(1).v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_neu__stimcodingdczprevresp;
            allresults(1).criterionshift    = results.repetition_neutral;
            allresults(1).subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');
            
            allresults(2).z_prevresp        = results.z_1_alt__stimcodingdczprevresp - results.z_2_alt__stimcodingdczprevresp;
            allresults(2).v_prevresp        = results.dc_1_alt__stimcodingdczprevresp - results.dc_2_alt__stimcodingdczprevresp;
            allresults(2).criterionshift    = results.repetition_alternating;
            allresults(2).subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');
            
            allresults(3).z_prevresp        = results.z_1_rep__stimcodingdczprevresp - results.z_2_rep__stimcodingdczprevresp;
            allresults(3).v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_rep__stimcodingdczprevresp;
            allresults(3).criterionshift    = results.repetition_repetitive;
            allresults(3).subjnr            = results.subjnr;
            alltitles{3}                    = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');
            
        case {'Bharath_fMRI', 'Anke_MEG', 'Anke_merged'}
            
            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY
            
            % use the stimcoding difference only from alternating
            allresults(1).z_prevresp        = results.z_1_0_50_0__stimcodingdczprevresp - results.z_2_0_50_0__stimcodingdczprevresp;
            allresults(1).v_prevresp        = results.dc_1_0_50_0__stimcodingdczprevresp - results.dc_2_0_50_0__stimcodingdczprevresp;
            allresults(1).criterionshift    = results.repetition_neutral;
            allresults(1).subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');
            
            allresults(2).z_prevresp        = results.z_1_0_20_0__stimcodingdczprevresp - results.z_2_0_20_0__stimcodingdczprevresp;
            allresults(2).v_prevresp        = results.dc_1_0_20_0__stimcodingdczprevresp - results.dc_2_0_20_0__stimcodingdczprevresp;
            allresults(2).criterionshift    = results.repetition_alternating;
            allresults(2).subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');
            
            allresults(3).z_prevresp        = results.z_1_0_80_0__stimcodingdczprevresp - results.z_2_0_80_0__stimcodingdczprevresp;
            allresults(3).v_prevresp        = results.dc_1_0_80_0__stimcodingdczprevresp - results.dc_2_0_80_0__stimcodingdczprevresp;
            allresults(3).criterionshift    = results.repetition_repetitive;
            allresults(3).subjnr            = results.subjnr;
            alltitles{3}                   = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');
            
        otherwise
            
            try
                % use the stimcoding difference
                results.z_prevresp = ...
                    results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;
            catch
                results.z_prevresp = ...
                    results.z_1_0__stimcodingdczprevresp - results.z_2_0__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1_0__stimcodingdczprevresp - results.dc_2_0__stimcodingdczprevresp;
            end
            
            results.criterionshift = results.repetition;
            
            % assign to structure
            allresults(1).z_prevresp = results.z_prevresp;
            allresults(1).v_prevresp = results.v_prevresp;
            allresults(1).criterionshift = results.criterionshift;

            alltitles{1} = datasetnames{d}{1}; % use only the dataset title
    end
    
    disp(datasets{d}); disp(numel(unique(results.subjnr)));
    close all;
    
    % PLOT
    subplot(4,4,1); hold on;
    [rho1, tt1] = plotScatter(allresults, 'z_prevresp', 0.585);
    xlabel('History bias in z', 'interpreter', 'tex', 'color', 'k');
    ylabel('P(repeat)');
    
    sp2 = subplot(4,4,2); hold on;
    [rho2, tt2] = plotScatter(allresults, 'v_prevresp', 0.05);
    xlabel('History bias in v', 'interpreter', 'tex', 'fontweight', 'normal', 'color', 'k');
    set(gca, 'yticklabel', []);
    
    % compute the difference in correlation
    [rho3, pval3] = corr(cat(1, allresults(:).v_prevresp), cat(1, allresults(:).z_prevresp), ...
        'rows', 'complete', 'type', 'pearson');
    if pval3 < 0.05,
        fprintf('warning %s: rho = %.3f, pval = %.3f \n', datasets{d}, rho3, pval3);
    end
    [rhodiff, ~, pval] = rddiffci(rho1,rho2,rho3,numel(~isnan( cat(1, allresults(:).criterionshift))), 0.05);
    offsetAxes; drawnow;
    
    % move together
    sp2.Position(1) = sp2.Position(1) - 0.08;
    ss = suplabel(datasetnames{d}{1}, 't');
    set(ss, 'fontweight', 'normal');
    ss.FontWeight = 'normal';
    ss.Position(2) = ss.Position(2) - 0.007;
    
    %% add line between the two correlation coefficients
    txt = {sprintf('\\Deltar_{%d} = %.3f, p = %.3f', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3, rhodiff, pval)};
    if pval < 0.001,
        txt = {sprintf('\\Deltar_{%d} = %.3f, p < 0.001', length(find(~isnan(cat(1, allresults(:).criterionshift) )))-3,  rhodiff)};
    end
    title(txt, 'fontweight', 'bold', 'fontsize', 6, 'horizontalalignment', 'left');
    
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1c_HDDM_modelfree_stimcoding_d%d.pdf', d));
end

end

function [rho, tt] = plotScatter(allresults, fld, legendWhere);

% overall correlation
x = cat(1, allresults(:).(fld));
y = cat(1, allresults(:).criterionshift);

% show line
[rho, pval] = corr(x,y, 'type', 'pearson', 'rows', 'complete');

% CORRELATION LINE SEPARATELY FOR EACH DATASET?
p = polyfit(x, ...
    y, 1);
xrange = linspace(min(x), max(x), 100);
yrange = polyval(p, xrange);
l = plot(xrange, yrange);
l.Color = 'k';
l.LineWidth = 0.5;
if pval < 0.05,
    l.LineStyle = '-';
else
    l.LineStyle = ':';
end

axis square;

% show lines to indicate origin
xlims = [min(x) max(x)];
ylims = [min(y) max(y)];
plot([0 0], ylims, 'color', [0.5 0.5 0.5], 'linewidth', 0.2);
plot(xlims, [0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 0.2); % if p(repeat), 0.5

% color in different grouos

% use separate colors for the different transition probabilities, from
% Thomas, blue green pink
% transitioncolors = [[8/256 141/256 165/256]; [165/256 8/256 141/256]; [141/256 165/256 8/256]];
transitioncolors = [[0.5 0.5 0.5]; [8/256 141/256 165/256]; [141/256 165/256 8/256]; [165/256 8/256 141/256]];
meancolors = [0 0 0; 0 0 1; 0 1 0];

markers = {'o', '^', 'v'}; %also indicate with different markers

for a = length(allresults):-1:1, % neutral last
    
    % PLOT ALL DATAPOINTS IN SPECIFIC COLOR
    s  = scatter(allresults(a).(fld), allresults(a).criterionshift,  10, ...
        markers{a}, 'LineWidth', 0.001, ...
        'markeredgecolor', 'w', 'markerfacecolor', transitioncolors(a, :));
end

for a = length(allresults):-1:1, % neutral last
    % also add the group mean
    p = ploterr(nanmean(allresults(a).(fld)), nanmean(allresults(a).criterionshift), 2*nanstd(allresults(a).(fld)) ./ sqrt(length(allresults(a).(fld))), ...
        2*nanstd(allresults(a).criterionshift) ./ sqrt(length(allresults(a).criterionshift)), '.', 'abshhxy', 0);
    set(p(1), 'markersize', 0.1, 'color', meancolors(a, :)); % tiny marker
    set(p(2), 'color', meancolors(a, :), 'linewidth', 0.5);
    set(p(3), 'color', meancolors(a, :), 'linewidth', 0.5);
end

axis tight; offsetAxes;

% PRINT THE CORRELATION COEFFICIENT
txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2, rho) sprintf('p = %.3f', pval)};
if pval < 0.001,
    txt = {sprintf('r_{%d} = %.3f', length(find(~isnan(y)))-2,rho) sprintf('p < 0.001')};
end
tt = text(min(get(gca, 'xlim')) + legendWhere*(range(get(gca, 'xlim'))), ...
    min(get(gca, 'ylim')) + 0.8*(range(get(gca, 'ylim'))), ...
    txt, 'fontsize', 5);
set(gca, 'color', 'none');

end
