function e1_serialBias_SfN_DIC()

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames


set(groot, 'defaultaxesfontsize', 5, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');
usr = getenv('USER');

types = {'stimcoding'};
for s = 1:length(types),
    plots = {'neutral', 'biased'};
    
    for p = 1:length(plots),
        
        % neutral vs biased plots
        switch p
            case 1
                
                switch usr
                    case 'anne' % local
                        datasets = {'RT_RDK', 'Anke_2afc_neutral', 'NatComm', 'projects/0/neurodec/Data/MEG-PL'};
                    case 'aeurai' % lisa/cartesius
                        datasets = {'RT_RDK', 'Anke_neutral', 'NatComm', 'MEG'};
                end
                datasetnames = { {'2AFC, RT', 'n = 22'}, {'2AFC, Braun et al. 2017 bioRxiv', 'n = 22'}, ...
                    {'2IFC, Urai et al. 2017 NatComm', 'n = 27'}, {'2IFC, replication', 'n = 61'}};
                
            case 2
                switch usr
                    case 'anne' % local
                        datasets = {'Anke_2afc_alternating', 'Anke_2afc_neutral', 'Anke_2afc_repetitive'};
                    case 'aeurai' % lisa/cartesius
                        datasets = {'Anke_alternating', 'Anke_neutral', 'Anke_repetitive', 'Anke_serial'};
                end
                datasetnames = {'2AFC alternating', '2AFC neutral', '2AFC repetitive', '2AFC all'};
                datasetnames =  {{'2AFC, Braun et al. 2017 bioRxiv', 'Alternating'}, ...
                    {'2AFC, Braun et al. 2017 bioRxiv', 'Neutral'}, ...w
                    {'2AFC, Braun et al. 2017 bioRxiv', 'Repetitive'}, ...
                    {'2AFC, Braun et al. 2017 bioRxiv', 'All'}};
        end
        
        % ============================================ %
        % DIC COMPARISON BETWEEN DC, Z AND BOTH
        % ============================================ %
        
        % 1. STIMCODING, only prevresp
        close all;
        mdls = {'dc_prevresp_prevstim', 'z_prevresp_prevstim', ...
            'dc_z_prevresp_prevstim', 'nohist'};
        for d = 1:length(datasets),
            subplot(4, 4, d);
            getPlotDIC(mdls, types{s}, d, 1);
            title(datasetnames{d});
            set(gca, 'xtick', 1:3, 'xticklabel', {'dc', 'z', 'both'});
            xlabel({'Modulation by', 'previous response & stimulus'});
        end
        
        % can you please add the full diffusion equation on top of the regression models? So akin to the eq we put into JW???s figure,
        % only here we would also have to add z. For specifics, take a look at the first few equations in Bogasz??? paper...
        
        drawnow; tightfig;
        print(gcf, '-depsc', sprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_%s_prevresp_prevstim.eps', plots{p}, types{s}));
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/figure1b_HDDM_DIC_%s_%s_prevresp_prevstim.pdf', plots{p}, types{s}));
        
        % % ============================================ %
        % DIC COMPARISON BETWEEN DC, Z AND BOTH
        % ============================================ %
        
        % build up models
        % 'regress_nohist' % should go last, will be baseline
        % 'regress_dc_prevresp'
        % 'regress_dc_prevresp_prevstim'
        % regress_dc_prevresp_prevstim_vasessions # add learning effects
        % then add modulation
        % then also add multiple responses into the past
        
        close all; nrsubpl = length(datasets);
        mdls = {'dc_prevresp', ...
            'dc_prevresp_prevstim',  ...
            'dc_prevresp_prevstim_prevrt', ...
            'dc_prevresp_prevstim_prevrt_prevpupil', ...
            'dc_prevresp_prevstim_vasessions',  ...
            'dc_prevresp_prevstim_vasessions_prevrt', ...
            'dc_prevresp_prevstim_vasessions_prevrt_prevpupil', ...
            'nohist'};
        
        for d = 1:length(datasets),
            subplot(nrsubpl, nrsubpl, d);
            getPlotDIC(mdls, types{s}, d, 1);
            set(gca, 'xtick', 1:length(mdls)-1);
            title(datasetnames{d});
        end
        
        switch types{s}
            case 'regress'
                
                subplot(nrsubpl, nrsubpl, nrsubpl+1);
                text(0, -0.2, {'Regression models', ...
                    '[1] v ~ 1 + stimulus + prevresp', ...
                    '[2] v ~ 1 + stimulus + prevresp + prevstim', ...
                    '[3] v ~ 1 + stimulus + prevresp*prevrt + prevstim*prevrt', ...
                    '[4] v ~ 1 + stimulus*session + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil', ...
                    '[5] v ~ 1 + stimulus*session + prevresp + prevstim', ...
                    '    a ~ 1 + session', ...
                    '[6] v ~ 1 + session*stimulus + prevresp*prevrt + prevstim*prevrt', ...
                    '     a ~ 1 + session', ...
                    '[7] v ~ 1 + session*stimulus + prevresp*prevrt + prevstim*prevrt + prevresp*prevpupil + prevstim*prevpupil', ...
                    '     a ~ 1 + session', ...
                    }, 'fontsize', 6); axis off;
        end
        
        tightfig;
        print(gcf, '-depsc', sprintf('~/Data/serialHDDM/suppfigure1b_HDDM_DIC_allmodels_%s_%s.eps', plots{p}, types{s}));
         print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/suppfigure1b_HDDM_DIC_allmodels_%s_%s.pdf', plots{p}, types{s}));
    end
end
end

% ============================================ %
% DIC COMPARISON BETWEEN DC, Z AND BOTH
% ============================================ %

function getPlotDIC(mdls, s, d, plotBest)

global datasets
axis square;

mdldic = nan(1, length(mdls));
for m = 1:length(mdls),
    if ~exist(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
            datasets{d}, s, mdls{m}), 'file'),
        continue;
    end
    
    load(sprintf('~/Data/%s/HDDM/summary/%s_%s_all.mat', ...
        datasets{d}, s, mdls{m}));
    
    if (isnan(dic.full) || isempty(dic.full)) && ~all(isnan(dic.chains)),
        dic.full = nanmean(dic.chains);
    end
    mdldic(m) = dic.full;
end

if isnan(mdldic(end)), assert(1==0); end

% everything relative to the full mdoel
mdldic = bsxfun(@minus, mdldic, mdldic(end));
mdldic = mdldic(1:end-1);

b = bar(1:length(mdldic), mdldic, 'facecolor', [0.8 0.8 0.8], 'barwidth', 0.5, 'BaseValue', 0);

% indicate which one is best
bestcolor = linspecer(3);
if plotBest,
    [~, idx] = min(mdldic);
    hold on;
    %  disp(idx);
    bar(idx, mdldic(idx), 'basevalue', 0, 'facecolor', [0.6 0.6 0.6], 'barwidth', 0.5, 'BaseValue', 0);
end

%# Add a text string above/below each bin
for i = 1:length(mdldic),
    if mdldic(i) < 0,
        text(i, mdldic(i) + 0.14*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 3, 'horizontalalignment', 'center');
    elseif mdldic(i) > 0,
        text(i, mdldic(i) + 0.14*range(get(gca, 'ylim')), ...
            num2str(round(mdldic(i))), ...
            'VerticalAlignment', 'top', 'FontSize', 3, 'horizontalalignment', 'center');
    end
end
xlim([0.5 length(mdldic)+0.5]);
ylabel('\Delta DIC (from nohist)', 'interpreter', 'tex');
offsetAxes; box off;
axis square;

end
