function e6_serialBias_SfN_modelFree_CRF_PPC
% ========================================== %
% conditional response functions from White & Poldrack
% run on simulated rather than real data
% ========================================== %

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames

qntls = [.2, .4, .6, .8, .95]; % White & Poldrack
qntls = [.1, .3, .5, .7, .9, 1]; % Leite & Ratcliff

for d = 3; % 1:length(datasets),
    
    % plot
    close all;
    subplot(441); hold on;
    % colors = flipud(cbrewer('div', 'PuOr', 5)); % repeaters in purple, alternators in orange
    plot(qntls, 0.5*ones(size(qntls)), 'color',  'k');
    
    % redo this for each simulation
    models = {'stimcoding_z_prevresp', 'stimcoding_dc_prevresp', 'stimcoding_nohist'};
    colors = [141 165 8;  8 141 165; 150 150 150] ./ 256;

    for m = 3:length(models),
        
        % load simulated data - make sure this has all the info we need
        alldata    = readtable(sprintf('~/Data/HDDM/%s/%s/ppc_data.csv', datasets{d}, models{m}));
        
        % use the simulations rather than the subjects' actual responses
        alldata.rt = alldata.rt_sampled; 
        alldata.response = alldata.response_sampled;
        
        % define repeaters and alternators based on dc
        dat = readtable(sprintf('~/Data/HDDM/summary/%s/allindividualresults.csv', datasets{d}));
        dat = dat(dat.session == 0, :);
        
        % recode into repeat and alternate
        alldata.repeat = zeros(size(alldata.response));
        alldata.repeat(alldata.response == (alldata.prevresp > 0)) = 1;
        
        % recode into biased and unbiased choices 
        alldata.biased = alldata.repeat;
        altIdx = ismember(alldata.subj_idx, unique(dat.subjnr(dat.repetition < 0.5)));
        alldata.biased(altIdx) = double(~(alldata.biased(altIdx))); % flip
        
        % divide RT into quantiles for each subject
        discretizeRTs = @(x) {discretize(x, quantile(x, [0, qntls]))};
        rtbins = splitapply(discretizeRTs, alldata.rt, findgroups(alldata.subj_idx));
        alldata.rtbins = cat(1, rtbins{:});
        % scatter(1:length(alldata.rt), alldata.rt, 10, alldata.rtbins, 'filled');
        
        % get RT quantiles for choices that are in line with or against the bias
        [gr, sjidx, rtbins] = findgroups(alldata.subj_idx, alldata.rtbins);
        cpres        = array2table([sjidx, rtbins], 'variablenames', {'subj_idx', 'rtbin'});
        cpres.choice = splitapply(@nanmean, alldata.biased, gr); % choice proportion
        
        % make into a subjects by rtbin matrix
        mat = unstack(cpres, 'choice', 'rtbin');
        mat = mat{:, 2:end}; % remove the last one, only has some weird tail
        
        % biased choice proportion
        plot(qntls, nanmean(mat, 1), 'color', colors(m, :), 'linewidth', 2);
    end
    
    axis tight; box off;
    set(gca, 'xtick', qntls);
    axis square;  offsetAxes;
    xlabel('RT (quantiles)'); ylabel('Fraction biased choices)');
    title(datasetnames{d}{1});
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/CRF_PPC_d%d.pdf', d));
    
end

end
