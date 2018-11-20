function compare_svgroups

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all;
global mypath datasets datasetnames 

colors = cbrewer('seq', 'Blues', 5);

% STARTING POINT POSTERIORS
close; cnt = 1;

for gr = [0 1 2],
    for d = 1:length(datasets),

        nohist = readtable(sprintf('%s/%s/stimcoding_nohist_svgroup/group_traces.csv', mypath, datasets{d}));
        withhist = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp_svgroup/group_traces.csv', mypath, datasets{d}));

        subplot(3,6,cnt); hold on; cnt = cnt + 1;

        ax1 = violinPlot_distribution(d, nohist.(['sv_' num2str(gr) '_']), [0 0 0]);
        ax2 = violinPlot_distribution(d, withhist.(['sv_' num2str(gr) '_']), colors(end, :));
        
        % add something to indicate significance
        difference = nohist.(['sv_' num2str(gr) '_']) - withhist.(['sv_' num2str(gr) '_']);

        p_neg = mean(difference > 0);
        p_pos = mean(difference < 0);
        mysigstar(gca, mean(get(gca, 'xlim')), ...
            min(get(gca, 'ylim')) + 0.1*range(get(gca, 'ylim')), min([p_neg p_pos]));
        title({datasetnames{d}{1}, datasetnames{d}{2}, sprintf('quantile %d', gr)});
    end
    end

    l = legend([ax1 ax2], {'no history', 'z and vbias'});
    l.Position(2) = l.Position(2) - 0.25;
    l.Box = 'off';
    suplabel('Drift-rate variability', 'y');

    disp('done')
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/svgroup_comparison.pdf', d));


% % SANITY CHECK
% close; cnt = 1;

% for gr = [0 1 2],
%     for d = 1:length(datasets),

%         %nohist = readtable(sprintf('%s/%s/stimcoding_nohist_svgroup/group_traces.csv', mypath, datasets{d}));
%         withhist = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp_svgroup/group_traces.csv', mypath, datasets{d}));

%         subplot(3,6,cnt); hold on; cnt = cnt + 1;

%         assert(1==0)
%         %ax1 = violinPlot_distribution(d, nohist.(['sv_' num2str(gr) '_']), [0 0 0]);
%         ax2 = violinPlot_distribution(d, withhist.(['sv_' num2str(gr) '_']), colors(end, :));
        
%         % % add something to indicate significance
%         % difference = nohist.(['sv_' num2str(gr) '_']) - withhist.(['sv_' num2str(gr) '_']);

%         % p_neg = mean(difference > 0);
%         % p_pos = mean(difference < 0);
%         % mysigstar(gca, mean(get(gca, 'xlim')), ...
%         %     min(get(gca, 'ylim')) + 0.1*range(get(gca, 'ylim')), min([p_neg p_pos]));
%         title({datasetnames{d}{1}, datasetnames{d}{2}, sprintf('quantile %d', gr)});
%     end
%     end

%     % l = legend([ax1 ax2], {'no history', 'z and vbias'});
%     % l.Position(2) = l.Position(2) - 0.25;
%     % l.Box = 'off';
%     suplabel('|History shift in v_{bias}|', 'y');

%     disp('done')
%     print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/svgroup_comparison_historybias.pdf', d));
% end
