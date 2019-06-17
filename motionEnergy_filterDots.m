function motionEnergy_filterDots(subjects, sessions, runs)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

% run with: sj = 1:33, sess = 1:3, run = 1:6
%
% requires code from
% Urai AE & Wimmer K (2016). Spatiotemporal motion energy filtering: 
% a Matlab implementation. Zenodo, http://doi.org/10.5281/zenodo.594505.

addpath('~/code/motionEnergy');
addpath(genpath('~/code/MEG'));
addpath(genpath('~/code/Tools'));
addpath('~/Documents/fieldtrip');
ft_defaults;

% get the dot coordinates for motion energy
% determine the path to the data
datapath = '/home/abraun/meg_data'; % home directory in home024/btalluri
savepath = '~/Data/motionEnergyAnkeMEG';

if nargin < 1,
    
    %% CONCATENATE INTO ONE HUGE TABLE
    files   = dir(sprintf('%s/*file*.mat', savepath));
    for f = 1:length(files),
        tmp = load(sprintf('%s/%s', files(f).folder, files(f).name));
        tmpdat.motionenergy{f} = tmp.motionenergy.trial;
        tmpdat.time{f} = tmp.motionenergy.time;
        tmpdat.behav{f} = tmp.behav;
    end
    
    data.motionenergy   = cat(1, tmpdat.motionenergy{:}); % flip, so that it positively scales with stimulus
    data.behavior       = cat(1, tmpdat.behav{:});
    assert(size(data.motionenergy, 1) == size(data.behavior, 1));
    data.timeaxis       = cat(1, tmpdat.time{:});
    data.timeaxis       = nanmean(data.timeaxis, 1);

    if ~istable(data.behavior), data.behavior = struct2table(data.behavior); end
    savefast(sprintf('%s/motionEnergyData_AnkeMEG.mat', savepath), 'data');
    return;
end

if ischar(subjects), subjects = str2double(subjects); end
if ischar(sessions), sessions = str2double(sessions); end
if ischar(runs), runs = str2double(runs); end

for sj = subjects,
    for sess = sessions,
        for run = runs,
            
            % stim is the dot coordinates and has the following dimensions: runs x trials x frames x 2 x dots
            % for each run, number of: runs = 1, trials = 40, frames = 45, dots = 117
            filename = dir(sprintf('%s/P%02d/Behav/MEG_P%d_s%d_r%d_*.mat', ...
                datapath,sj,sj,sess,run));
            % assert(length(filename) == 1, 'could not find a unique file');
            
            for fileidx = 1:length(filename),
                
                resultsFile = sprintf('%s/P%02d_session%d_run%d_file%d.mat', ...
                    savepath,sj,sess,run,fileidx);
                
                % don't rerun
                if exist(resultsFile, 'file'),
                    fprintf('%s already exists, skipping\n', savepath,resultsFile);
                   % return;
                end
                
                load(sprintf('%s/P%02d/Behav/%s', datapath,sj, filename(fileidx).name),...
                    'stim','dots','results','window','setup');
                disp(filename(fileidx).name);
                
                % avoid window UI error
                display.dist        = window.dist;
                display.res         = window.res;
                display.width       = window.width;
                display.frameRate   = window.frameRate;
                display.center      = window.center;
                display.ppd         = deg2pix(display, 1);
                coord               = stim;
                
                % make sure the framerate matches the nr of dot frames, 750 ms dots
                assert(abs(0.75 - size(coord, 3) / display.frameRate) < 0.02, ...
                    'frameRate in dots and behav does not match!');
                
                % ==================================================================
                % BUILD FILTERS
                % ==================================================================
                
                % temporal range of the filter
                cfg            = [];
                cfg.frameRate  = display.frameRate;
                cfg.ppd        = display.ppd;
                
                % k = 60, from Kiani et al. 2008
                cfg.k = 60;
                
                % adjust spatial filters to match the speed in the dots
                effectiveSpeed = pix2deg(display, dots.speed) ./ dots.nvar;
                
                % Kiani et al. 2008 has a speed of 2.84 deg/s and used sigma_c and sigma_g
                % as 0.35 (not explicitly mentioned in their paper). To optimally scale the
                % filters for the speed in our dots, multiply the spatial parameters
                % see email exchange with Klaus Wimmer
                cfg.sigma_c = 0.35 * (effectiveSpeed / 2.84);
                cfg.sigma_g = 0.35 * (effectiveSpeed / 2.84);
                
                % equations exactly as in Kiani et al. 2008
                [f1, f2] = makeSpatialFilters(cfg);
                [g1, g2] = makeTemporalFilters(cfg);
                
                % ==================================================================
                % get behavioural info
                % ==================================================================
                
                % start loop over trials
                clear data;
                ntrials         = size(coord, 2);
                data.trial      = nan(ntrials, size(coord, 3));
                
                % ==================================================================
                % FILTER EACH TRIAL
                % ==================================================================
                
                for t = 1:ntrials,
                    
                    % GET DOT COORDINATES FOR THIS TRIAL
                    thiscoord = squeeze(coord(1, t, :, :, :));
                    
                    % check that we have no NaNs left
                    assert(~any(isnan(thiscoord(:))));
                    
                    % change coordinates to movie
                    % always make the dots go up, that way will get the opponent
                    % energy for up-down
                    thisstim = coord2stim(display, thiscoord, 90);
                    
                    % filters
                    motionenergy = applyFilters(thisstim, f1, f2, g1, g2);
                    
                    % also save all the energy summed across space
                    data.trial(t, :) = squeeze(sum(sum(motionenergy)));
                    
                end
                
                % time axis
                data.time = 0 : 1/display.frameRate : (size(motionenergy, 3) - 1)/display.frameRate;
                assert(numel(data.time) == size(thisstim, 3));
                
                % ==================================================================
                % SANITY CHECK
                % ==================================================================
                
                % for the last one, check that this worked;
                singletrial_motionenergy = nanmean(data.trial, 2);
                roc = rocAnalysis(singletrial_motionenergy(dots.direction == 270), ...
                    singletrial_motionenergy(dots.direction == 90), 0, 0);
                if roc.i < 0.9,
                    warning('motion energy does not separate stimulus types');
                end
                
                % save this block
                motionenergy = data;
                
                %% ADD BEHAVIORAL VARIABLES
                behav = array2table([sign(dots.direction - 100)' dots.coherence' ...
                    sign(results.response - 100)' ...
                    results.correct' results.RT', transpose(1:length(results.response))], ...
                    'variablenames', ...
                    {'stimulus', 'coherence', 'response', 'correct', 'RT', ...
                    'trialnum'});
                behav.prevresp      = circshift(behav.response, 1);
                behav.prevstim      = circshift(behav.stimulus, 1);
                behav.prevrt        = circshift(behav.stimulus, 1);
                behav.prevcorrect   = circshift(behav.correct, 1);
                behav.prev2resp     = circshift(behav.response, 2);
                behav.prev2stim     = circshift(behav.stimulus, 2);
                behav.prev2correct  = circshift(behav.correct, 2);
                behav.subj_idx      = sj *ones(size(behav.response));
                behav.session       = sess *ones(size(behav.response));
                behav.block         = run *ones(size(behav.response));
                behav.transitionprob = setup.transition_probability * ones(size(behav.response));
                behav               = table2struct(behav);
                
                save(resultsFile, 'motionenergy', 'behav');
                disp(resultsFile);
                
                %                 behav.motionenergy = motionenergy.trial;
                %                 writetable(behav, sprintf('%s/P%02d_session%d_run%d.csv', ...
                %                     savepath,sj,sess,run));
            end
        end
    end
end

end