function multiplicative_vbias_psychfuncs_ppc

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath colors

% redo this for each simulation
models = {'stimcoding_dc_z_prevresp', 'stimcoding_dc_z_prevresp_multiplicative', 'data'};
% BLACK FOR DATA, BLUE FOR DC_Z, CYAN FOR MULTIPLICATIVE
cyan = cbrewer('seq', 'GnBu', 10);
thesecolors = {mean(colors([1 2], :)), cyan(5, :)};

for d = [5]
    
    close all; subplot(441); hold on;
    for m = 1:length(models),
        
        switch models{m}
            case 'data'
                filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
                alldata  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
            otherwise
                if ~exist(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}), 'file'),
                    continue;
                else
                    fprintf('%s/summary/%s/%s_ppc_data.csv \n', mypath, datasets{d}, models{m});
                end
                % load simulated data - make sure this has all the info we need
                alldata             = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}));
                alldata.rt          = abs(alldata.rt_sampled);
                alldata.response    = alldata.response_sampled;
        end
        
        % in the NCOMMS data, group the 3 easiest coherence levels
        if d == 5,
            alldata.coherence = alldata.coherence * 100;
            alldata.coherence(alldata.coherence < 5) = 2.5;
        end
        alldata.stimulus = alldata.coherence .* sign(alldata.stimulus-0.01);
        
        % make sure to use absolute RTs!
        alldata             = sortrows(alldata, {'subj_idx'});
        alldata.rt          = abs(alldata.rt);
        
        % recode into repeat and alternate for the model
        alldata.repeat      = zeros(size(alldata.response));
        alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
        
        % for each observers, compute their bias
        [gr, sjs] = findgroups(alldata.subj_idx);
        sjrep = splitapply(@nanmean, alldata.repeat, gr);
        sjrep = sjs(sjrep < 0.5);
        altIdx = ismember(alldata.subj_idx, sjrep);
        
        % recode into biased and unbiased choices
        alldata.pref = alldata.prevresp;
        alldata.pref(altIdx) = -alldata.pref(altIdx);
        alldata.pref(alldata.pref == -1) = 0;
        
        alldata.biased = alldata.repeat;
        alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
        
        %         % get RT quantiles for choices that are in line with or against the bias
        %         [gr, sjidx, cohs, biased] = findgroups(alldata.subj_idx, alldata.stimulus, alldata.pref);
        %         cpres               = array2table([sjidx, cohs biased], 'variablenames', {'subj_idx', 'coh', 'biased'});
        %         cpres.choice        = splitapply(@nanmean, alldata.response, gr); % choice proportion
        %
        %         % make into a psychometric functions
        %         mat         = unstack(cpres, 'choice', 'coh');
        %         mat_pref    = mat{mat.biased == 1, 3:end};
        %         mat_unpref  = mat{mat.biased == 0, 3:end};
        %
        %         %% FIRST, PANEL WITH PREFERRED AND UNPREFERRED PSYCHOMETRIC FUNCTIONS
        %         switch models{m}
        %             case 'data'
        %                 colors2 = cbrewer('qual', 'Set1', 9);
        %                 colors2 = colors2([9 4], :);
        %
        %                 close all; subplot(441); hold on;
        %                 errorbar(1:length(unique(cohs)), nanmean(mat_pref), ...
        %                     nanstd(mat_pref) ./ sqrt(size(mat_pref, 1)), ...
        %                     's-', 'capsize', 0, 'markerfacecolor', 'w', 'color', colors2(1, :), 'markersize', 2);
        %                 errorbar(1:length(unique(cohs)), nanmean(mat_unpref), ...
        %                     nanstd(mat_unpref) ./ sqrt(size(mat_pref, 1)), ...
        %                     'd-', 'capsize', 0, 'markerfacecolor', 'w', 'color', colors2(2, :), 'markersize', 2);
        %                 xlabel('Sensory evidence (%)');
        %                 ylabel('Choice fraction');
        %                 set(gca, 'xtick', 1:length(unique(cohs)), 'xticklabel', unique(cohs));
        %                 box off; % axis square;
        %                 title(datasetnames{d});
        %                 offsetAxes; tightfig;
        %                 print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/psychfunc_biased_%d.pdf', d));
        %                 close all; subplot(441); hold on;
        %         end
        %
        %         %% COMPUTE BIAS AS A FUNCTION OF ABSOLUTE COHERENCE
        %         mat_biased = mat_pref - mat_unpref;
        %         switch d
        %             case 4
        %                 mat_biased = mean(cat(3, mat_biased(:, 5:-1:1), mat_biased(:, 5:end)), 3);
        %             case 5
        %                 mat_biased =  mean(cat(3, mat_biased(:, 5:-1:1), mat_biased(:, 6:end)), 3);
        %         end
        
        
        clear cpres;        
        % get RT quantiles for choices that are in line with or against the bias
        [gr, sjidx, cohs] = findgroups(alldata.subj_idx, abs(alldata.stimulus));
        cpres               = array2table([sjidx, cohs], 'variablenames', {'subj_idx', 'coh'});
        cpres.choice        = splitapply(@nanmean, alldata.biased, gr); % choice proportion
        
        % make into a psychometric functions
        mat         = unstack(cpres, 'choice', 'coh');
        mat_biased = mat{:, 2:end};

        switch models{m}
            case 'data'

                errorbar(1:length(unique(abs(alldata.stimulus))), nanmean(mat_biased), ...
                    nanstd(mat_biased) ./ sqrt(size(mat_biased, 1)), ...
                    'o-', 'capsize', 0, 'markerfacecolor', 'k', 'color', 'k', 'markersize', 3);
            otherwise
                %  plot(1:length(unique(abs(alldata.stimulus))), nanmean(mat_biased), '-', 'color', thesecolors{m});
                boundedline(1:length(unique(abs(alldata.stimulus))), nanmean(mat_biased), ...
                    nanstd(mat_biased) ./ sqrt(size(mat_biased, 1)), 'cmap', thesecolors{m}, 'alpha');
                
        end
        
        mat_biased_save{m} = mat_biased;
    end
    
    % compare the two model predictions
    %[h, pval] = ttest(mat_biased_save{3}, mat_biased_save{1}, 'alpha', 0.05/length(unique(abs(alldata.stimulus))));
    %[h, pval] = ttest(mat_biased_save{3}, mat_biased_save{1}, 'alpha', 0.05/length(unique(abs(alldata.stimulus))));
    
    
    xlabel('Evidence strength (%)');
    ylabel('P(bias)');
    set(gca, 'xtick', 1:length(unique(abs(alldata.stimulus))), 'xticklabel', unique(abs(alldata.stimulus)));

    box off; % axis square;
    title(datasetnames{d});
    
    axis tight;
    set(gca,  'ylim', [0.5 max(get(gca, 'ylim'))]);
    offsetAxes; tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/psychfunc_biased_ppc_%d.pdf', d));
    
end
