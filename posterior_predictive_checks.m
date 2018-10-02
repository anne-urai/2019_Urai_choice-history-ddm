function posterior_predictive_checks

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

addpath(genpath('~/code/Tools'));
warning off; close all; clear;
global datasets datasetnames mypath

% ========================================== %
% MODULATION OF SERIAL CHOICE BIAS
% ========================================== %

% datasets = {'JW_PNAS', 'JW_yesno', 'Murphy', 'Anke_MEG_neutral', 'NatComm', 'MEG'};
plotWhich = 'stimcoding'; % {'error', 'biased', 'stimcoding'};
choiceCat = {{'left', 'right'}, {'down', 'up'}, {'weaker', 'stronger'}, {'weaker', 'stronger'}, {'no', 'yes'}, {'no','yes'}};

for d = 1:length(datasets),
    close all;
    
    if ~exist(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, 'stimcoding_nohist'), 'file'),
        fprintf('cannot find %s/stimcoding_nohist/ppc_data.csv \n', datasets{d});
        continue;
    else
        disp(datasets{d});
    end
    
    % get traces for the model with pupil and rt modulation
    ppc = readtable(sprintf('%s/summary/%s/%s_ppc_data.csv', mypath, datasets{d}, 'stimcoding_nohist'));
    ppc.correct                        = (ppc.stimulus == ppc.response);
    ppc.repeat                         = zeros(size(ppc.response));
    ppc.repeat(ppc.response == (ppc.prevresp > 0)) = 1;
    
    % for each observers, compute their bias
    [gr, sjs] = findgroups(ppc.subj_idx);
    sjrep = splitapply(@nanmean, ppc.repeat, gr);
    sjrep = sjs(sjrep < 0.5);
    
    % recode real data into biased vs unbiased
    ppc.biased                         = ppc.repeat;
    altIdx                             = ismember(ppc.subj_idx, sjrep);
    ppc.biased(altIdx) = double(~(ppc.biased(altIdx))); % flip
    
    switch plotWhich
        case 'error'
            ppc.biased = ppc.correct;
            ppc.response = ppc.biased;
        case 'stimcoding'
            ppc.biased = (ppc.response > 0);
    end
    
    % unbiased RTs negative
    ppc.rt(ppc.biased == 1)           = abs(ppc.rt(ppc.biased == 1));
    ppc.rt(ppc.biased == 0)           = -abs(ppc.rt(ppc.biased == 0));
    
    % SAME FOR THE SIMULATED DATA
    ppc.correct_sampled = (ppc.stimulus == ppc.response_sampled);
    
    % recode into repeat and alternate for the model
    ppc.repeat_sampled = zeros(size(ppc.response_sampled));
    ppc.repeat_sampled(ppc.response_sampled == (ppc.prevresp > 0)) = 1;
    
    % recode into biased and unbiased choices
    ppc.biased_sampled = ppc.repeat_sampled;
    altIdx = ismember(ppc.subj_idx, sjrep);
    ppc.biased_sampled(altIdx) = double(~(ppc.biased_sampled(altIdx))); % flip
    
    switch plotWhich
        case 'error'
            ppc.biased_sampled = ppc.correct_sampled;
        case 'stimcoding'
            ppc.biased_sampled = (ppc.response_sampled > 0);
    end
    
    % define the sampled RT also by the sampled correctness!
    ppc.modelcorrect                   = (ppc.response_sampled == ppc.stimulus);
    ppc.rt_sampled(ppc.biased_sampled == 1)   = abs(ppc.rt_sampled(ppc.biased_sampled == 1));
    ppc.rt_sampled(ppc.biased_sampled == 0)   = -abs(ppc.rt_sampled(ppc.biased_sampled == 0));
    ppc = ppc(:, {'rt', 'rt_sampled', 'stimulus', 'response'}); % save some memory
    
    % determine the colors
    switch plotWhich
        case 'error'
            bestcolor = cbrewer('qual', 'Dark2', 5);
            bestcolor = bestcolor([2 1], :);
            fitcolor = cbrewer('qual', 'Set2', 5);
            fitcolor = fitcolor([2 1], :);
        case 'biased'
            bestcolor = cbrewer('div', 'PiYG', 6);
            bestcolor = bestcolor([1 end], :);
            fitcolor = [0 0 0];
        case 'stimcoding'
            bestcolor = cbrewer('qual', 'Dark2', 5);
            bestcolor = bestcolor([3 5], :);
            fitcolor = cbrewer('qual', 'Set2', 5);
            fitcolor = fitcolor([3 5], :);
    end
    
    switch plotWhich
        case {'error', 'biased'};
            ppc.stimulus = ones(size(ppc.stimulus));
    end
    
    ix = unique(ppc.stimulus);
    rx = unique(ppc.response);
    for i = 1:length(ix),
        sph{i} = subplot(4,10,i);
        hold on;
        for r = 1:length(rx),
            histogram_smooth(abs(ppc.rt(ppc.stimulus == ix(i) & ppc.response == rx(r))), ...
                abs(ppc.rt_sampled(ppc.stimulus==ix(i) & ppc.response == rx(r))), ...
                bestcolor(r, :), bestcolor(r, :), fitcolor(r, :));
        end
        
        axis tight; % axis square;
        offsetAxes_y;
        maxRT = ceil(max(abs(ppc.rt)));
        if maxRT == 5, maxRT = 4; end
        if maxRT < 3, maxRT = 3; end

        % if d > 3,maxRT = 3; end
        xlim([0 maxRT]); set(gca, 'xtick', [0 maxRT], 'xminortick', 'on');
        %  ylabel('Probability');
        switch plotWhich
            case 'stimcoding'
        title({'Stimulus', capitalize(choiceCat{d}{i})}, 'color', bestcolor(i, :), 'fontweight', 'normal');
        end
        set(gca, 'yticklabel', []);
        set(gca, 'xcolor', 'k', 'ycolor', 'k');
        
    end
    
    try
        % move together
        sph{2}.Position(1) = sph{2}.Position(1) - 0.01;
    end
    
    % xlabel('RT (s)');
        [ss, h1] = suplabel('RT (s)', 'x');
        ss.Position(2) = ss.Position(2) + 0.04;
        h1.Color = 'k';
    [ss, h1] = suplabel('Probability', 'y');
    ss.Position(1) = ss.Position(1) + 0.06;
    h1.Color = 'k';
    try
        set(sph{2}, 'ylim', get(sph{1}, 'ylim'));
    end
    
    % legend for choices!
    switch plotWhich
        case 'stimcoding'
            ylims = get(gca, 'ylim');
            text(maxRT*0.7, max(ylims)*0.7, 'Choice', 'fontsize', 6);
            text(maxRT*0.7, max(ylims)*0.6, sprintf('"%s"', capitalize(choiceCat{d}{1})), 'color', bestcolor(1, :), 'fontsize', 6);
            text(maxRT*0.7, max(ylims)*0.5, sprintf('"%s"', capitalize(choiceCat{d}{2})), 'color', bestcolor(2, :), 'fontsize', 6);
            
            %set(gcf, 'color', 'none');
            set(gca, 'xcolor', 'k', 'ycolor', 'k');
            %[ss, h1] =  suplabel(cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2}), 't');
            %ss.Position(2) = ss.Position(2) + 0.04;
    end
    
    tightfig;
    switch plotWhich
        case 'error'
            print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d.pdf', d));
        case 'biased'
            print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d_biased.pdf', d));
        case 'stimcoding'
            print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/PPC_d%d_stimcode.pdf', d));
            % export_fig(sprintf('~/Data/serialHDDM/PPC_d%d_stimcode.eps', d));
    end
    
    
end

end

function h = histogram_smooth(x1, x2, color1, color2, fitcolor)

% % manually count so i can plot myself
% [n, edges] = histcounts(x1, -3:0.1:3, 'normalization', 'pdf');
%
% posidx = find(edges > 0); posidx(posidx > length(n)) = [];
% negidx = find(edges < 0);
%
% % plot as stairs??
% %bar(edges(posidx), n(posidx), 'edgecolor', 'none', 'facecolor', color1, 'barwidth', 1);
% %bar(edges(negidx), n(negidx), 'edgecolor', 'none', 'facecolor', color2, 'barwidth', 1);
%
% % [n, edges] =

% first the fit - make sure this is not normalized to 1!
%[f,xi] = ksdensity(x2);
%h = plot(xi, f, 'color', fitcolor, 'linewidth', 0.75);

% put the real number of trials on the y-axis

[n1, edges1] = histcounts(x1, -5:0.05:5);
[n2, edges2] = histcounts(x2, -5:0.05:5); % much smaller steps, smoother

% correctionRatio
%n2 = n2*10;
stairs(edges1(1:end-1), n1, 'color', fitcolor, 'linewidth', 1);
plot(edges2(1:end-1), n2, 'color', color1, 'linewidth', 0.75);

% histogram(x2, -3:0.01:3, 'displaystyle', 'stairs', ...
%     'edgecolor', fitcolor, 'linewidth', 0.75);
%
% % put the real number of trials on the y-axis
% histogram(x1, -3:0.1:3, 'displaystyle', 'stairs', ...
%     'edgecolor', color1, 'linewidth', 0.75);

% remove white box in the pdf
set(gca, 'color', 'none');

end

function offsetAxes_y()

if ~exist('ax', 'var'), ax = gca;
end
if ~exist('offset', 'var'), offset = 4;
end

% ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/offset;
ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/offset;

% this will keep the changes constant even when resizing axes
addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax));

end

function resetVertex ( ax )
% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
ax.XRuler.Axle.VertexData(1,2) = max(get(ax, 'Xtick'));
end