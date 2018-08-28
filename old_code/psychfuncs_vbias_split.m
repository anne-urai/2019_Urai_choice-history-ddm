function psychfuncs_vbias_split

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath colors

% redo this for each simulation
models = {'stimcoding_nohist', 'stimcoding_z_prevresp',  ...
    'stimcoding_dc_prevresp', 'stimcoding_dc_z_prevresp' ...
    'data'};
thesecolors = {[0.5 0.5 0.5],  colors(1, :), ...
    colors(2, :), mean(colors([1 2], :)), [0 0 0]};

for RTs = [0 1 2];
    for d = [4 5],
        close all;
        
        % plot
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
                    alldata    = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, models{m}));
                    alldata    = sortrows(alldata, {'subj_idx'});
            end
            
            if m < length(models),
                % use the simulations rather than the subjects' actual responses
                alldata.rt          = abs(alldata.rt_sampled);
                alldata.response    = alldata.response_sampled;
            end
            
            % make sure to use absolute RTs!
            alldata.rt = abs(alldata.rt);
            alldata.correct = ((alldata.stimulus > 0) == alldata.response);
            assert(nanmean(alldata.correct) > 0.5);
            
            %% DEPENDING ON RT FLAG, USE ALL/FAST/SLOW TRIALS
            switch RTs
                case 0
                    % do nothing, keep all
                case 1
                    % keep fast RTs
                    for sj = unique(alldata.subj_idx)',
                        for c = unique(abs(alldata.coherence))',
                            trls = find(alldata.subj_idx == sj & abs(alldata.coherence) == c);
                            medRT = nanmedian(alldata.rt(trls));
                            % remove slow RTs
                            alldata.rt(alldata.subj_idx == sj & abs(alldata.coherence) == c & ...
                                alldata.rt > medRT) = NaN;
                        end
                    end
                case 2
                    % keep slow RTs
                    for sj = unique(alldata.subj_idx)',
                        for c = unique(abs(alldata.coherence))',
                            trls = find(alldata.subj_idx == sj & abs(alldata.coherence) == c);
                            medRT = nanmedian(alldata.rt(trls));
                            % remove fast RTs
                            alldata.rt(alldata.subj_idx == sj & abs(alldata.coherence) == c & ...
                                alldata.rt < medRT) = NaN;
                        end
                    end
            end
            alldata(isnan(alldata.rt), :) = [];
            
            % recode into repeat and alternate for the model
            alldata.repeat = zeros(size(alldata.response));
            alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
            
            % for each observers, compute their bias
            [gr, sjs] = findgroups(alldata.subj_idx);
            sjrep = splitapply(@nanmean, alldata.repeat, gr);
            sjrep = sjs(sjrep < 0.5);
            
            % recode into biased and unbiased choices
            alldata.biased = double(alldata.prevresp > 0);
            altIdx = ismember(alldata.subj_idx, sjrep);
            alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
            % biased [0,1] now codes for repeat [0,1] for repeaters and switch
            % [0,1] for alternaters
            
            if numel(unique(alldata.stimulus)) > 2,
                alldata.stimulus = sign(alldata.stimulus);
            elseif isequal(unique(alldata.stimulus), [0 1]'),
                alldata.stimulus = sign(alldata.stimulus - 0.1);
            end
            
            % from fraction to percent
            if all(abs(alldata.coherence) < 1),
                alldata.coherence = alldata.coherence * 100;
            end
            
            alldata.coherence       = [alldata.coherence .* alldata.stimulus];
            [gr, sjidx, biased, coh] = findgroups(alldata.subj_idx, alldata.biased, alldata.coherence);
            psychfunc               = array2table([sjidx, biased, coh], 'variablenames', {'subj_idx', 'pref', 'coh'});
            psychfunc.resp          = splitapply(@nanmean, alldata.response, gr); % choice proportion
            
            % make into a subjects x coherence x biased matrix
            mat = unstack(psychfunc, 'resp', 'subj_idx');
            subplot(3,3,m); hold on;
            
            % FIT 4-PARAMETER PSYCHOMETRIC FUNCTIONS
            mat0 = mat{mat.pref == 0, 3:end};
            mat1 = mat{mat.pref == 1, 3:end};
            for sj = 1:length(unique(sjidx)),
                [params1(1, sj) params1(2, sj) params1(3, sj) params1(4, sj)] = ...
                    fitLogistic(unique(coh), mat0(:, sj));
                [params2(1, sj) params2(2, sj) params2(3, sj) params2(4, sj)] = ...
                    fitLogistic(unique(coh), mat1(:, sj));
            end
            
            % compare parameters between these two types of trials
            paramdiff = nanmean(params1 - params2, 2);
            [h, p] = ttest(params1', params2');
            
            % read out the curve and plot
            logisticY = @(x,p) p(3)+(1-p(3)-p(4)) * (1./(1+exp(- ( p(1) + p(2).*x ))));
            xx = linspace(min(coh), max(coh), 100);
            y1 = logisticY(xx, nanmean(params1, 2));
            y2 = logisticY(xx, nanmean(params2, 2));
            plot(xx, y1, 'color', thesecolors{m}, 'linestyle', '-');
            plot(xx, y2, 'color', thesecolors{m}, 'linestyle', ':');
            
            % PLOT DATAPOINTS ON TOP
            h = ploterr(unique(coh), nanmean(mat{mat.pref == 0, 3:end}, 2), [], ...
                1.96 *  nanstd(mat{mat.pref == 0, 3:end}, 0, 2) ./ sqrt(numel(unique(sjidx))), 'k', 'abshhxy', 0);
            set(h(1), 'color', thesecolors{m}, 'marker', 'o', ...
                'markerfacecolor', 'w', 'markeredgecolor', thesecolors{m}, 'linewidth', 0.5, 'markersize', 2, ...
                'linestyle', 'none');
            set(h(2), 'linewidth', 0.5, 'color', thesecolors{m});
            
            h = ploterr(unique(coh),nanmean(mat{mat.pref == 1, 3:end}, 2), [], ...
                1.96 *  nanstd(mat{mat.pref == 1, 3:end}, 0, 2) ./ sqrt(numel(unique(sjidx))), 'k', 'abshhxy', 0);
            set(h(1), 'color', thesecolors{m}, 'marker', 'd', ...
                'markerfacecolor', 'w', 'markeredgecolor', thesecolors{m}, 'linewidth', 0.5, 'markersize', 2, ...
                'linestyle', 'none');
            set(h(2), 'linewidth', 0.5, 'color', thesecolors{m});
            
            % LAYOUT
            title(regexprep(models{m}, 'stimcoding_', ''), 'interpreter', 'none');
            axis tight;
            ylim([0 1]); box off;
            set(gca, 'xtick', sort([unique(coh(abs(coh) > 0.5*median(abs(coh)))); 0]));
            axis square;  offsetAxes;
            xlabel('Coherence (%)');
            set(gca, 'xcolor', 'k', 'ycolor', 'k');
            ylabel('P(response == 1)');
            
            for pix = 1:4,
                if p(pix) < 1e-3
                    txt = '***';
                elseif p(pix) < 1e-2
                    txt = '**';
                elseif p(pix) < 0.05
                    txt = '*';
                elseif ~isnan(p(pix)),
                    txt = 'n.s.';
                end
                pval{pix} = txt;
            end
            
            % ADD STATS TEXT
            disptxt = {sprintf('\\Delta\\mu %.2f, %s', paramdiff(1), pval{1}), ...
                sprintf('\\Delta\\sigma %.2f, %s', paramdiff(2), pval{2}), ...
                sprintf('\\Delta\\gamma %.2f, %s', paramdiff(3), pval{3}), ...
                sprintf('\\Delta\\lambda %.2f, %s', paramdiff(4), pval{4})};%, ...
            text(nanmean(abs(coh)), 0.3, disptxt, 'interpreter', 'tex');
        end
        
        suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
        tightfig;
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/psychfuncs_biased_d%d_RTs%d.pdf', d, RTs));
    end
end

end
