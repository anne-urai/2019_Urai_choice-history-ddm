function multiplicative_vbias_residuals_scatter

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

close all; clc;
addpath(genpath('~/code/Tools'));

global mypath datasets datasetnames

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

whichMeasures = {'mse', 'corr', 'resid'};
for w = 1:length(whichMeasures),
    
    models = {'stimcoding_dc_z_prevresp', 'stimcoding_dc_z_prevresp_multiplicative'};
    for d = [4 5] % only NatComm and Anke_MEG_neutral, with varying coherence level
        figure;
        
        for m = 1:length(models),
            
            disp(datasets{d});
            
            % load simulated data - make sure this has all the info we need
            alldata    = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}));
            alldata    = sortrows(alldata, {'subj_idx'});
            
            % recode into repeat and alternate for the model
            alldata.biased = ((alldata.response > 0) == (alldata.prevresp > 0));
            alldata.biased_model = ((alldata.response_sampled > 0) == (alldata.prevresp > 0));
            
            % for each observers, compute their bias
            [gr, sjs] = findgroups(alldata.subj_idx);
            sjrep = splitapply(@nanmean, alldata.biased, gr);
            sjrep = sjs(sjrep < 0.5);
            
            % recode into biased and unbiased choices
            altIdx                      = ismember(alldata.subj_idx, sjrep);
            alldata.biased(altIdx)      = double(~(alldata.biased(altIdx))); % flip
            alldata.biased_model(altIdx) = double(~(alldata.biased_model(altIdx))); % flip
            
            % FOR BOTH REAL AND SIMULATED CHOICES, SEE BIAS PROPORTION AS A
            % FUNCTION OF COHERENCE
            [gr, sj, coh]   = findgroups(alldata.subj_idx, alldata.coherence);
            tab = array2table([sj, coh], 'variablenames', {'subj_idx', 'coherence'});
            tab_data = tab;
            tab_model = tab;
            tab_data.bias   = splitapply(@nanmean, alldata.biased, gr);
            tab_model.bias  = splitapply(@nanmean, alldata.biased_model, gr);
            mat_data        = unstack(tab_data, 'bias', 'coherence');
            mat_data = mat_data{:, 2:end};
            mat_model       = unstack(tab_model, 'bias', 'coherence');
            mat_model = mat_model{:, 2:end};
            
            cohtab = repmat(unique(coh)', length(unique(sj)), 1);
            
            colors  = cbrewer('seq', 'PuBuGn', numel(unique(coh)) + 5);
            colors  = colors([3:end-4 end], :);
            colormap(colors);
            %
            %         subplot(441);
            %         scatter(mat_data(:), mat_model(:), 10, cohtab(:));
            %         xlabel('P(bias) observed');
            %         ylabel('P(bias) predicted');
            %         title(datasetnames{d});
            %         axis square; axis tight;  axisEqual; offsetAxes;
            %
            %         tightfig;
            %         print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence_residuals_scatter_d%d_m%d.pdf', d, m));
            %
            %% mean squared error
            
            switch whichMeasures{w}
                case 'mse'
                    mse = nanmean((mat_model-mat_data).^2);
                case 'corr'
                    mse = diag(corr(mat_model, mat_data, 'rows', 'pairwise'));
                case 'resid'
                    mse = nanmean(mat_data-mat_model);
            end
            
            subplot(441); hold on;
            if m == 1,
                plot(1:length(unique(coh)), mse, '-', 'color', [0 0 0]);
                s1 =  scatter(1:length(unique(coh)), mse, 20, unique(coh), 'filled');
                s1.MarkerEdgeColor = 'k';
                %  scatter(1:length(unique(coh)), mse, 25, [0.3 0.3 0.3]);
                
            elseif m == 2,
                plot(1:length(unique(coh)), mse, '-', 'color', [0 0 0 ]);
                s2 = scatter(1:length(unique(coh)), mse, 20, unique(coh), 's', 'filled');
                s2.MarkerEdgeColor = 'k';
                
                %  scatter(1:length(unique(coh)), mse, 25, [0 0 0]);
            end
            
        end
        
        switch whichMeasures{w}
            case 'corr'
                ylabel({'Correlation observed', 'vs. predicted P(bias)'});
            case 'mse'
                ylabel({'MSE observed', 'vs. predicted P(bias)'});
            case 'resid'
                ylabel({'Predicted -', 'observed P(bias)'});
        end
        set(gca, 'xtick', 1:length(unique(coh)), 'xticklabel', unique(coh)*100);
        
        axis tight; axis square;
        xlim([min([0.5 get(gca, 'xlim')]) max(get(gca, 'xlim'))]);
        
        offsetAxes;
        % legend([s1 s2], {'Single v_{bias}', 'Difficulty-dependent v_{bias}'});
        
        switch d
            case 4
                xlabel('Coherence (%)');
            case 5
                xlabel('\Delta coherence (%)');
                set(gca, 'xticklabelrotation', -30);
        end
        % title(datasetnames{d});
        
        tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/vbias_per_coherence_residuals_mse_d%d_m%d_w%d.pdf', d, m, w));
        
    end
end


