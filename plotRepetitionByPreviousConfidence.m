function plotRepetitionByPreviousConfidence
% Anne Urai, 2018
% anne.urai@gmail.com

clear all; clc; 
nbins = 2;

%% determine how the figures will look
set(groot, 'defaultaxesfontsize', 6, 'defaultaxestitlefontsizemultiplier', 1.1, ...
    'defaultaxeslabelfontsizemultiplier', 1.1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.05], 'defaultaxestickdir', 'out', 'DefaultAxesTickDirMode', 'manual', ...
    'defaultfigurecolormap', [1 1 1], 'defaultTextInterpreter','tex');

usr = getenv('USER');
switch usr
    case 'anne'
        mypath = '~/Data/HDDM';
    case 'aeurai'
        mypath  = '/nfs/aeurai/HDDM';
end

datasets = {'Murphy', 'JW_PNAS', 'JW_yesno', 'NatComm', 'MEG', 'Anke_2afc_sequential', 'Anke_MEG'};
datasetnames = { {'Murphy 2014'},  {'de Gee 2014'}, {'de Gee 2017'}, ...
    {'Urai 2017'}, {'Urai in prep'},  {'Braun 2018'}, {'Braun in prep'}, {'Talluri in prep'}};

addpath(genpath('~/Desktop/code/gramm/'));
warning off; close all;
cmap = viridis(256);
colormap(cmap);

whichVars = {'prev_rt', 'prev_difficulty'};
for v = 1:length(whichVars),
    
    dattotal = {}; close all;
    for d = 1:length(datasets),
        
        close all;
        csv_file = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
        dat = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, csv_file.name));
        
        % only use neutral
        if ismember('transitionprob', dat.Properties.VariableNames)
            dat = dat(dat.transitionprob == 0.5, :);
        end
        
        % for previous coherence, recode
        if strcmp(whichVars{v}, 'prev_difficulty'),
            if ismember('coherence', dat.Properties.VariableNames) || length(unique(dat.stimulus)) > 2,
                if ~ismember('coherence', dat.Properties.VariableNames),
                    dat.coherence = abs(dat.stimulus);
                end
                dat.prevcoh = circshift(dat.coherence, 1);
                dat.prevrt  = dat.prevcoh;
            else
                continue;
            end
            lb1 = 'hard'; lb2 = 'easy';
            axlb = 'Previous difficulty';
        else
            lb1 = 'fast'; lb2 = 'slow';
            axlb = 'Previous RT';
        end
        
        %% COMPUTE REPETITION PROBABILITY AS A FUNCTION OF PREVIOUS RT
        % TODO: WITHIN EACH PARTICIPANT!
        dat              = sortrows(dat, 'subj_idx');
        zscore_bin       = @(x) {discretize(x, [-inf quantile(x, nbins), inf])};
        tmp              = splitapply(zscore_bin, dat.prevrt, findgroups(dat.subj_idx));
        dat.prev_rt_bin  = vertcat(tmp{:});
        
        dat.repeat       = ((dat.response > 0) == (dat.prevresp > 0));
        dat.prev_correct = (dat.prevstim == dat.prevresp);
        
        [gr, prevRT, prevCorrect, sj] = findgroups(dat.prev_rt_bin, dat.prev_correct, dat.subj_idx);
        dat2 = array2table([prevRT, prevCorrect, sj], 'variablenames', {'prev_rt', 'prev_correct', 'subj_idx'});
        dat2.repeat = splitapply(@nanmean, dat.repeat, gr);
        
        %% SPLIT INTO REPEATERS AND ALTERNATORS
        summarydat       = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
        dat2.subgroup    = repmat({'Repeaters'}, height(dat2), 1);
        alternators      = unique(summarydat.subjnr(summarydat.repetition < 0.5));
        dat2.subgroup(ismember(dat2.subj_idx, alternators)) = {'Alternators'};
        dat2.dataset     = repmat(datasetnames{d}, height(dat2), 1);
        dat2.datasetnr   = repmat(d, height(dat2), 1);
        
        %% ALSO INCLUDE A NORMALIZED RT
        dat2             = sortrows(dat2, 'subj_idx');
        normalize        = @(x) {x - nanmean(x) + 0.5};
        repeat_norm      = splitapply(normalize, dat2.repeat, findgroups(dat2.subj_idx));
        repeat_norm_mat  = vertcat(repeat_norm{:});
        dat2.repeat_norm = repeat_norm_mat;
        dattotal{end+1}  = dat2;
    end
    
    %% ADD A GRAND AVERAGE SUBJECT
    dattotal        = vertcat(dattotal{:});
    [gr, prevRT, prevCorrect] = findgroups(dattotal.prev_rt, dattotal.prev_correct);
    dat3            = array2table([prevRT, prevCorrect], 'variablenames', {'prev_rt', 'prev_correct'});
    dat3.repeat     = splitapply(@nanmean, dattotal.repeat_norm, gr);
    dat3.repeat_norm     = splitapply(@nanmean, dattotal.repeat_norm, gr);
    dat3.subj_idx   = repmat(0, height(dat3), 1);
    dat3.dataset    = repmat({'Average'}, height(dat3), 1);
    dat3.datasetnr  = repmat(0, height(dat3), 1);
    dat3.subgroup   = repmat({'All'}, height(dat3), 1);
    dattotal        = vertcat(dattotal, dat3);  % merge
    
    %% USE GRAMM TO PLOT
    g = gramm('x', dattotal.prev_rt,'y', dattotal.repeat, 'color', dattotal.prev_correct);
    g.facet_wrap(dattotal.dataset, 'ncols',4);
    % g.geom_hline('yintercept', 0.5, 'style', '-', 'linewidth', 0.5);
    % Plot linear fits of the data with associated confidence intervals
    g.stat_summary('setylim', 1);
    % Set appropriate names for legends
    g.set_names('column', '', 'x', axlb, 'y',...
        'P(repeat)', 'color','Previous correct');
    % Do the actual drawing
    lb = repmat({' '}, 1, nbins-1);
    g.axe_property('PlotBoxAspectRatio', [1 1 1], 'ylim', [0.45 0.65], ...
        'ytick', 0.45:0.05:0.65, ...
        'xlim', [1 nbins+1], 'xtick', 1:nbins+1, 'xticklabel', {lb1, lb{:}, lb2}); % axis square
    g.draw();
    
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperOrientation','landscape');
    g.export('file_name', sprintf('~/Data/serialHDDM/%s_repetition.pdf', whichVars{v}), 'file_type', 'pdf');
end
end