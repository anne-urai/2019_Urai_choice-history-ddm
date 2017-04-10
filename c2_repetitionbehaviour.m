% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
warning off;

% ============================================ %
% TWO DIFFERENT DATASETS
% ============================================ %

usr = getenv('USER');
switch usr
    case 'anne' % local
        datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL'};
    case 'aeurai' % lisa/cartesius
        datasets = {'RT_RDK', 'MEG-PL'};
end

set(groot, 'defaultaxesfontsize', 10, 'defaultaxestitlefontsizemultiplier', 1, ...
    'defaultaxestitlefontweight', 'bold', ...
    'defaultfigurerenderermode', 'manual', 'defaultfigurerenderer', 'painters');

for d = 2:length(datasets),
    results = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', datasets{d}));
    
    % ============================================ %
    % compute repetition parameters from separate HDDM models
    % ============================================ %
    
    results.dc_prevresp__stimcoding_prevresp_dc = ...
        results.dc_1__stimcoding_prevresp_dc - results.dc_2__stimcoding_prevresp_dc;
    
    results.z_prevresp__stimcoding_prevresp_z = ...
        results.z_2__stimcoding_prevresp_z - results.z_1__stimcoding_prevresp_z;
    
    results.dc_prevresp__stimcoding_prevresp_dc_z = ...
        results.dc_1__stimcoding_prevresp_dc_z - results.dc_2__stimcoding_prevresp_dc_z;
    
    results.z_prevresp__stimcoding_prevresp_dc_z = ...
        results.z_2__stimcoding_prevresp_dc_z - results.z_1__stimcoding_prevresp_dc_z;
    
    % center around zero
    results.p_repeat = results.repetition - 0.5;
    
    % ============================================ %
    % RENAME PARAMETERS
    % ============================================ %
    
    results.Properties.VariableNames{'dc_prevresp__stimcoding_prevresp_dc'}     = 'dc_seq_stimcoding';
    results.Properties.VariableNames{'v_prevresp__regress_dc_prevresp'}         = 'dc_seq_regress';
    results.Properties.VariableNames{'z_prevresp__stimcoding_prevresp_z'}       = 'z_seq_stimcoding';
    results.Properties.VariableNames{'z_prevresp__regress_z_prevresp'}          = 'z_seq_regress';
    results.Properties.VariableNames{'dc_prevresp__stimcoding_prevresp_dc_z'}   = 'dc_seq_stimcoding_joint';
    results.Properties.VariableNames{'v_prevresp__regress_dc_z_prevresp'}       = 'dc_seq_regress_joint';
    results.Properties.VariableNames{'z_prevresp__stimcoding_prevresp_dc_z'}    = 'z_seq_stimcoding_joint';
    results.Properties.VariableNames{'z_prevresp__regress_dc_z_prevresp'}       = 'z_seq_regress_joint';
    
    % ============================================ %
    % SEPARATE OR JOINT FIT
    % ============================================ %
    
    close;
    corrplot(results, {'dc_seq_stimcoding', ...
        'z_seq_stimcoding'}, ...
        {'dc_seq_stimcoding_joint', ...
        'z_seq_stimcoding_joint'});
    suplabel('Separate fits', 'x');
    suplabel('Joint fit', 'y');
    suplabel('Stimcoding', 't');
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/separate_vs_jointfit_stimcoding.pdf', datasets{d}));
    
    close;
    corrplot(results, {'dc_seq_regress', ...
        'z_seq_regress'}, ...
        {'dc_seq_regress_joint', ...
        'z_seq_regress_joint'});
    suplabel('Separate fits', 'x');
    suplabel('Joint fit', 'y');
    suplabel('Regression models', 't');
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/separate_vs_jointfit_regression.pdf', datasets{d}));
    
    % ============================================ %
    % STIMCODING VS REGRESSION MODELS
    % ============================================ %
    
    close;
    corrplot(results, {'dc_seq_stimcoding_joint', ...
        'z_seq_stimcoding_joint'}, ...
        {'dc_seq_regress_joint', ...
        'z_seq_regress_joint'});
    suplabel('Stimcoding', 'x');
    suplabel('Regression', 'y');
    suplabel('Joint fits', 't');
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/stimcoding_vs_regression.pdf', datasets{d}));
    
    close; corrplot(results, {'dc_seq_stimcoding', ...
        'z_seq_stimcoding'}, ...
        {'dc_seq_regress', ...
        'z_seq_regress'});
    suplabel('Stimcoding', 'x');
    suplabel('Regression', 'y');
    suplabel('Separate fits', 't');
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/stimcoding_vs_regression_separate.pdf', datasets{d}));
    
    % ============================================ %
    % CORRELATION WITH MODELFREE BEHAVIOUR
    % ============================================ %
    
    % use only those that are simultaneously fitted
    close;
    corrplot(results, {'repetition',  'criterionshift', ...
        'dc_seq_regress_joint', ...
        'z_seq_regress_joint'});
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/corr_with_modelfree.pdf', datasets{d}));
    
    % ============================================ %
    % CORRELATION WITH MODELFREE BEHAVIOUR
    % ============================================ %
    
    % use only those that are simultaneously fitted
    close;
    corrplot(results, {'p_repeat',  'criterionshift', ...
        'dc_seq_regress_joint', ...
        'z_seq_regress_joint'});
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/corr_with_modelfree.pdf', datasets{d}));
    
    % ============================================ %
    % STABILITY OF SERIAL CHOICE BIAS
    % ============================================ %
    
    subjects    = unique(results.subjnr(results.session == 2));
    
    results.repeat_s1 = nan(size(results.p_repeat));
    results.repeat_s1(results.session == 0 & ismember(results.subjnr, subjects)) = ...
        results.p_repeat(results.session == 1 & ismember(results.subjnr, subjects));
    results.repeat_s2 = nan(size(results.p_repeat));
    results.repeat_s2(results.session == 0 & ismember(results.subjnr, subjects)) = ...
        results.p_repeat(results.session == 2 & ismember(results.subjnr, subjects));
    
    results.criterionshift_s1 = nan(size(results.p_repeat));
    results.criterionshift_s1(results.session == 0 & ismember(results.subjnr, subjects)) = ...
        results.criterionshift(results.session == 1 & ismember(results.subjnr, subjects));
    results.criterionshift_s2 = nan(size(results.p_repeat));
    results.criterionshift_s2(results.session == 0 & ismember(results.subjnr, subjects)) = ...
        results.criterionshift(results.session == 2 & ismember(results.subjnr, subjects));
    
    close;
    corrplot(results, {'repeat_s1', ...
        'criterionshift_s1'}, ...
        {'repeat_s2', ...
        'criterionshift_s2'});
    suplabel('Session 1', 'x');
    suplabel('Session 2', 'y');
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/stability.pdf', datasets{d}));
    
    % ============================================ %
    % use separate HDDM fits
    % ============================================ %
    
    results2 	= readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults_separatesessions.csv', datasets{d}));
    results2.dc_seq__stimcoding_prevresp_dc_z = results2.dc_1__stimcoding_prevresp_dc_z - results2.dc_2__stimcoding_prevresp_dc_z;
    results2.z_seq__stimcoding_prevresp_dc_z = results2.z_2__stimcoding_prevresp_dc_z - results2.z_1__stimcoding_prevresp_dc_z;
    
    vars        = results2.Properties.VariableNames';
    subjects    = unique(results2.subjnr(results2.session == 2));
    
    clf; cnt = 1;
    for v = [31 35 36 39 40],
        subplot(3,3,cnt); cnt = cnt + 1;
        plot(results2.(vars{v})(results2.session == 1 & ismember(results2.subjnr, subjects)), ...
            results2.(vars{v})(results2.session == 2 & ismember(results2.subjnr, subjects)), '.');
        axis square; axisNotSoTight; box off;
        xlabel(vars{v}(1:end-26), 'interpreter', 'none');
        ylabel(vars{v}(1:end-26), 'interpreter', 'none');
        [rho, pval] = corr(results2.(vars{v})(results2.session == 1 & ismember(results2.subjnr, subjects)), ...
            results2.(vars{v})(results2.session == 2 & ismember(results2.subjnr, subjects)), 'type', 'spearman', 'rows', 'complete');
        if pval < 0.05, lsline; end
        title(sprintf('\\rho %.2f p %.3f', rho, pval));
    end
    suplabel('Session 1, MEG', 'x');
    suplabel('Session 2, MEG', 'y');
    
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/stability_HDDM.pdf', datasets{d}));
    
    % ============================================ %
    % DIC VALUES
    % ============================================ %
    
    models = {'regress_dc_prevresp', 'regress_z_prevresp', 'regress_dc_z_prevresp', ...
        'stimcoding_prevresp_dc', 'stimcoding_prevresp_z', 'stimcoding_prevresp_dc_z'};
    alldic = nan(30, length(models));
    for m = 1:length(models),
        load(sprintf('~/Data/%s/HDDM/summary/%s_all.mat', datasets{d}, models{m}));
        alldic(:, m) = (dic);
    end
    alldic = nanmean(alldic);
    
    close; subplot(221);
    bar(alldic(1:3), 'basevalue', nanmean(alldic(1:3)), 'facecolor', linspecer(1));
    box off;
    ylabel('DIC');
    set(gca, 'xticklabel', {'dc', 'z', 'both'});
    title('Regression models');
    set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
    axisNotSoTight;
    print(gcf, '-dpdf', sprintf('~/Data/%s/HDDM/summary/DICcomparison.pdf', datasets{d}));
    
    subplot(222);
    bar(alldic(4:6), 'basevalue', nanmean(alldic(4:6)), 'facecolor', linspecer(1));
    box off;
    ylabel('DIC');
    set(gca, 'xticklabel', {'dc', 'z', 'both'});
    title('Stimcoding');
    set(gca, 'YTickLabel', num2str(get(gca, 'YTick')'));
    axisNotSoTight;
    
    % ============================================ %
    % ADD RT AND PUPIL IN
    % ============================================ %
    
    % RT and pupil
    results.Properties.VariableNames{'v_prevpupil_prevresp__regress_dc_prevresp_prevpupil_prevrt'} ...
        = 'dc_pupil_seq_regress_joint';
    results.Properties.VariableNames{'v_prevrt_prevresp__regress_dc_prevresp_prevpupil_prevrt'} ...
        = 'dc_rt_seq_regress_joint';
    results.Properties.VariableNames{'z_prevpupil_prevresp__regress_z_prevresp_prevpupil_prevrt'} ...
        = 'z_pupil_seq_regress_joint';
    results.Properties.VariableNames{'z_prevrt_prevresp__regress_dc_z_prevresp_prevpupil_prevrt'} ...
        = 'z_rt_seq_regress_joint';
    
    results.dc_seq_regress   = results.v_prevresp__regress_dc_prevresp_prevpupil_prevrt;
    results.z_seq_regress    = results.z_prevresp__regress_z_prevresp_prevpupil_prevrt;
    
    close;
    corrplot(results, ...
        {'dc_seq_regress_joint'}, ...
        {'dc_rt_seq_regress_joint', 'dc_pupil_seq_regress_joint'});
    
    close;
    corrplot(results, ...
        {'z_seq_regress_joint'}, ...
        {'z_rt_seq_regress_joint', 'z_pupil_seq_regress_joint'});
    
end
