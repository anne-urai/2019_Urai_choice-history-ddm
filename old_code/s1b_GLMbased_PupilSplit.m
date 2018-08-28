function [] = s1b_GLMbased_PupilSplit()

addpath(genpath('~/Dropbox/code/RT_RDK'));
addpath('~/Documents/MATLAB/toolboxes/fieldtrip/');
ft_defaults;
clc; dbstop if error;
close all; clear all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.canonicalIRF      = true; % use canonical vs individual IRF params
cfg.IRFderivative     = true; % add derivative of IRF

cfg.runGLM            = false; % run the model (if false, just plotting)
cfg.timelock          = false; % add additional timelocking of the predicted data, takes more time

if ~exist('subjects', 'var'),
    % get all those subjects who have a folder
    cd ~/Data/RT_DDM/;
    s = dir('P*');
    s = {s([s(:).isdir]).name};
    for i = 1:length(s), subjects(i) = str2num(s{i}(2:3)); end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if cfg.timelock,
    stimfig = figure; respfig = figure; fbfig = figure;
end

if cfg.runGLM, % otherwise, only plotting
    
    for sj = (unique(cfg.subjects)),
        
        disp(sj);
        load(sprintf('~/Data/RT_DDM/P%02d_alleye_continuous.mat', sj));
        
        % get all the data
        fulltimecourse = data.trial{1}';
        trl = data.trialinfo;
        
        % remove trials without a response
        trl(find(isnan(trl(:,4))), :) = [];
        
        % remove RT outliers
        RT      = trl(:,6);
        trl(find(RT > 2.99), :) = [];
        trl(find(RT < 0.15), :) = [];
        
        % correct weirdness
        trl(find(trl(:,2)==27), 2) = 270;
        trl(find(trl(:,2)==9), 2) = 90;
        trl(find(trl(:,4)==27), 4) = 270;
        trl(find(trl(:,4)==9), 4) = 90;
        
        assert(all(trl(:,2)>50), 'wrongly coded stimulus');
        assert(all(trl(:,4)>50), 'wrongly coded stimulus');
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare regressors
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rcnt = 1;
        regressors(rcnt).name          = 'stim';
        regressors(rcnt).shape         = 'stick';
        regressors(rcnt).timecols      = 3;
        rcnt = rcnt + 1;
        
        regressors(rcnt).name          = 'stim-resp box';
        regressors(rcnt).shape         = 'box';
        regressors(rcnt).timecols      = [3 7];
        rcnt = rcnt + 1;
        
        regressors(rcnt).name          = 'resp';
        regressors(rcnt).shape         = 'stick';
        regressors(rcnt).timecols      = 7;
        rcnt = rcnt + 1;
        
        regressors(rcnt).name          = 'resp-fb ramp';
        regressors(rcnt).shape         = 'box';
        regressors(rcnt).timecols      = [7 9];
        rcnt = rcnt + 1;
        
        regressors(rcnt).name          = 'feedback';
        regressors(rcnt).shape         = 'stick';
        regressors(rcnt).timecols      = 9;
        
        %% BOOTSTRAP, MAXIMIZE THE DIFFERENCE BETWEEN TWO SETS OF BOXCARS
        
        nboot = 1000;
        
        % preallocate
        bootstrap.idx   = nan(nboot, length(trl));
        bootstrap.diff  = nan(nboot, 1);
        
        for boot = 1:nboot,
            
            % GET DESIGN MATRIX
            [designM, bootstrap.idx] = makeDesignMatrix(regressors, fulltimecourse, trl);
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % convolve with pupil IRF
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % choose whether to use canonical or not, and whether to add the
            % derivative or not
            [designM] = convolveIRF(designM, data.fsample, cfg.IRFderivative, cfg.canonicalIRF, sj);
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % run regression to get beta weights
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp('running GLM...');
            tic; [b, bint, ~, ~, stats] = regress(fulltimecourse, designM);  toc;
            
            % combine the betas using pythagoras
            if cfg.IRFderivative,
                nregressors = length(regressors);
                
                % check if the derivative is not bigger than the irf itself - this
                % will lead to strange results (talk to JW)
                for r = 1:nregressors,
                    if any(b(:, r) - b(:, r+nregressors)) > 0, warning('derivative is larger than IRF itself!'); end
                end
                
                % use pythagoras to take the length of b and b'
                b       = sqrt(b(:, 1:nregressors).^2 + b(:, nregressors+1:end).^2)';
                bint    = sqrt(bint(:, 1:nregressors).^2 + bint(:, nregressors+1:end).^2)';
                
            end
            
            % make sure this model makes any sense at all
            assert(stats(1) > 0, 'negative rsquared!');
            
            % get the measure we're interested in: the difference between
            % the high and low pupil split beta
            bootstrap.diff(boot) = b(2) - b(3);
            
        end
       
        % make an overview plot of this split
        figure; subplot(221);
        histogram(bootstrap.diff);
        
        % find the best split
        [maxdiff, best] = max(bootstrap.diff);
        bestidx = bootstrap.idx(best, :);
        
        % also do a timelock on those trials to see what the GLM is picking up
        
        subplot(2,2,2); hold on;
        for group = 0:1,
            
            % first, take subset of trials
            cfg         = [];
            cfg.trl     = trl(bestidx==group, :);
            data        = ft_redefinetrial(cfg, data);
            
            resplocked_data = shiftoffset_timelock( data, [], ...
                data.trialinfo(:,7)-data.trialinfo(:,3), 1, 3, data.fsample, 0);
            % plot
            boundedline(resplocked_data.time, [resplocked_data.avg resplocked_predicted.avg], ...
                permute([resplocked_data.var resplocked_predicted.var], [1 3 2]));
            
        end

    end
    
    save(sprintf('~/Data/RT_DDM/GrandAverage/pupilglm_%s_canonicalf%d_irfD%d.mat', cfg.whichrpe, cfg.canonicalIRF, cfg.IRFderivative), 'grandavg','cfg', 'regressors');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAND AVERAGE ACROSS THE GROUP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visualiseGrandAverage(cfg);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION THAT MAKES A DESIGN MATRIX
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [designM, idx] = makeDesignMatrix(regressors, fulltimecourse, trl)
% based on a set of regressors, output a design matrix
% EXTRA: RANDOMLY SPLIT UP THE STIM-RESP BOX INTO TWO SETS OF REGRESSORS

% determine how many regressors there will be,
% none of them are scaled by uncertainty here
nregressors = length(regressors);

% preallocate designM
designM     = zeros(length(fulltimecourse), nregressors);

% CREATE EACH REGRESSOR TIMEPOINT BY TIMEPOINT
colcnt = 1;
for r = 1:length(regressors),
    
    % see which type we have
    switch regressors(r).shape,
        case 'stick'
            designM(trl(:, regressors(r).timecols), colcnt) = 1;
        case 'box'
            
            % special case for this analysis
            switch regressors(r).name,
                case 'stim-resp box'
                    ntrials = length(trl);
                    p = randperm(ntrials);
                    idx = p > ntrials / 2; % binarize into 2 classes
                    
                    % make 2 sets of regressors here
                    for t = 1:size(trl, 1), 
                        if idx(t),
                            designM( trl(t, regressors(r).timecols(1)) : ...
                                trl(t, regressors(r).timecols(2)), colcnt) = 1;
                        else
                            designM( trl(t, regressors(r).timecols(1)) : ...
                                trl(t, regressors(r).timecols(2)), colcnt+1) = 1;
                        end
                    end
                    colcnt = colcnt + 1; % add extra col
                    
                otherwise
                    % normal box for all trials
                    for t = 1:size(trl, 1),  designM( trl(t, regressors(r).timecols(1)) : ...
                            trl(t, regressors(r).timecols(2)), colcnt) = 1;
                    end
            end
            
        case 'ramp'
            for t = 1:size(trl, 1),  designM( trl(t, regressors(r).timecols(1)) : ...
                    trl(t, regressors(r).timecols(2)), colcnt) = ...
                    linspace(0, 1, trl(t, regressors(r).timecols(2)) - trl(t, regressors(r).timecols(1)) + 1 );
            end
    end
    colcnt = colcnt + 1; % move one col to the right
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION THAT CONVOLVES DESIGNM WITH IRF
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [designM] = convolveIRF(designM, fsample, derivative, canonical, sj)
% makes a convolved version of the design matrix
% can use either canonical or individual parameters, and include the
% derivative of the IRF or not

if canonical,
    w = 10.1; tmax = 0.930; % from Hoeks and Levelt
else
    load('~/Data/UvA_pupil/GrandAverage/individualIRFs.mat');
    w    = grandavg.params(sj, 1);
    tmax = grandavg.params(sj, 2);
end

% Erlang gamma function from Hoeks and Levelt, Wierda
t = 0:1/(fsample):4;
irf = t.^w .* exp(-t.*w ./ tmax);

designM2 = nan(size(designM, 1) + length(irf) - 1, size(designM, 2));

% convolve the timepoints of each regressor
for r = 1:size(designM, 2),
    designM2(:, r) = conv(designM(:, r), irf, 'full');
end

% make sure we only take the part that matches with the original timecourse
designM2 = designM2(1:size(designM, 1), :);
assert(~any(isnan(designM2(:))));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add a set of regressors for the derivative of the irf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if derivative,
    
    irfD = diff(irf);
    designM3 = nan(size(designM, 1) + length(irfD) - 1, size(designM, 2));
    
    % convolve the timepoints of each regressor
    for r = 1:size(designM, 2),
        designM3(:, r) = conv(designM(:, r), irfD, 'full');
    end
    
    % make sure we only take the part that matches with the original timecourse
    designM3 = designM3(1:size(designM, 1), :);
    assert(~any(isnan(designM3(:))));
    
    designM = [ones(size(designM, 1), 1) designM2 designM3];
    
else
    % only use the normal IRF
    designM = [ones(size(designM, 1), 1) designM2];
end

if 0,
    % ZSCORE EACH REGRESSOR, SINCE THAT'S WHAT I DID WITH THE DATA
    % doesnt this introduce other artefacts again?
    for c = 2:size(designM, 2),
        designM(:, c) = zscore(designM(:, c));
    end
end

end
