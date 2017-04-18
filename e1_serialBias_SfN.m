function e1_serialBias_SfN()

addpath(genpath('~/code/Tools'));
warning off;
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'Anke_2afc_sequential'};
datasetnames = {'RT', '2IFC', 'Anke'};

set(groot, 'defaultaxesfontsize', 7, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

% 1. STIMCODING
clf;
mdls = {'stimcoding_dc_prevresp_prevstim', 'stimcoding_z_prevresp_prevstim', ...
    'stimcoding_dc_z_prevresp_prevstim'};
for d = 1:length(datasets),
    subplot(4,4,d); dic = nan(1, length(mdls));
    
    for m = 1:length(mdls),
        fileid = sprintf('~/Data/%s/HDDM/%s/DIC-combined.txt', ...
            datasets{d}, mdls{m});
        if exist(fileid, 'file'),
            txtfile = fopen(fileid);
            txtread = textscan(txtfile, 'Model %d: %d');
            dic(m)  = txtread{2};
            fclose(txtfile);
        end
        
        % in case the DIC could not be computed (custom link functions)
        if isnan(dic(m)),
            chains = 0:29;
            for c = 1:length(chains),
                fileid = sprintf('~/Data/%s/HDDM/%s/DIC-md%d.txt', ...
                    datasets{d}, mdls{m}, chains(c));
                if exist(fileid, 'file'),
                    txtfile = fopen(fileid);
                    txtread = textscan(txtfile, 'Model %d: %d');
                    dicchains(c)  = txtread{2};
                    fclose(txtfile);
                end
            end
            try
                dic(m) = mean(dicchains);
            end
        end
    end
    
    dic = dic - dic(3); % relative to full model
    try
        bar(dic(1:2), 'basevalue', dic(3), 'facecolor', linspecer(1));
    end
    box off;
    ylabel('\Delta DIC (from full model)');
    set(gca, 'xticklabel', {'dc', 'z'});
    title(datasetnames{d});
    % set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
end
suplabel('Stimcoding', 'x');

% 2. REGRESSION MODELS
mdls = {'regress_dc_prevresp_prevstim', 'regress_z_prevresp_prevstim', ...
    'regress_dc_z_prevresp_prevstim'};
for d = 1:length(datasets),
    subplot(4,4,d+8); dic = nan(1, length(mdls));
    
    for m = 1:length(mdls),
        fileid = sprintf('~/Data/%s/HDDM/%s/DIC-combined.txt', ...
            datasets{d}, mdls{m});
        if exist(fileid, 'file'),
            txtfile = fopen(fileid);
            txtread = textscan(txtfile, 'Model %d: %d');
            dic(m)  = txtread{2};
            fclose(txtfile);
        end
        
        % in case the DIC could not be computed (custom link functions)
        if dic(m) == 0,
            clear dicchains;
            chains = 0:29;
            for c = 1:length(chains),
                fileid = sprintf('~/Data/%s/HDDM/%s/DIC-md%d.txt', ...
                    datasets{d}, mdls{m}, chains(c));
                if exist(fileid, 'file'),
                    txtfile = fopen(fileid);
                    txtread = textscan(txtfile, 'Model %d: %d');
                    dicchains(c)  = txtread{2};
                    fclose(txtfile);
                end
            end
            try
                dic(m) = mean(dicchains);
            end
        end
    end
    
    dic = dic - dic(3); % relative to full model
    try
        bar(dic(1:2), 'basevalue', dic(3), 'facecolor', linspecer(1));
    end
    box off;
    ylabel('\Delta DIC (from full model)');
    set(gca, 'xticklabel', {'dc', 'z'});
    title(datasetnames{d});
    % set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
end
suplabel('Regression models', 'x');
print(gcf, '-dpdf', '~/Data/serialHDDM/DIC.pdf');

% ============================================ %
% SANITY CHECK; DOES DRIFT RATE CORRELATE WITH D'?
% ============================================ %

% ANKE'S DATA WITH DIFFERENT COHERENCE LEVELS

clf;
mdls = {'regress_dc_prevresp_prevstim_sessions'};
for d = 1:length(datasets),
    % load data
    csvfile = dir(sprintf('~/Data/%s/HDDM/*.csv', datasets{d}));
    alldata = readtable(sprintf('~/Data/%s/HDDM/%s', datasets{d}, csvfile.name));
    % compute a bunch of basic things from Matlab
    results     = b3b_behaviouralMetrics(alldata);
    
    % get point estimates from results
    hddmresults = readtable(sprintf('~/Data/%s/HDDM/%s/results-combined.csv', datasets{d}, mdls{1}));
    
    subplot(3,3,d);
    scatter(results.dprime(results.session == 0), hddmresults.dprime);
    
end
% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

