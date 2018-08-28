function [] = PupilAnalysisMain_Elisa()

addpath(genpath('~/Dropbox/code/Pupil_Learning'));
addpath('~/Documents/MATLAB/fieldtrip/');
ft_defaults;
clc; dbstop if error;

% get all those subjects who have a folder
cd ~/Data/RT_DDM/;
s = dir('P*');
s = {s(:).name};
for i = 1:length(s), subjects(i) = str2num(s{i}(2:3)); end

subjects = 12;

for sj = unique(subjects),
    tic;
    
    cd(sprintf('~/Data/RT_DDM/P%02d/', sj));
    
    if ~exist(sprintf('P%02d_alleye.mat', sj), 'file'),
        
        clear sessions;
        % check which sessions to use
        s = dir('S*');
        s = {s(:).name};
        for i = 1:length(s), sessions(i) = str2num(s{i}(2)); end
        
        for session = unique(sessions),
            
            cd(sprintf('~/Data/RT_DDM/P%02d/S%d', sj, session));
            
            % ==================================================================
            % LOAD IN SUBJECT SPECIFICS AND READ DATA
            % ==================================================================
            
            clear blocks;
            % check which blocks to use
            s = dir('P*.asc');
            s = {s(:).name};
            for i = 1:length(s), blocks(i) = sscanf(s{i},'P%*d_s%*d_%*[a-z]%d_*.asc'); end
            
            blocks = 1;
            
            for block = unique(blocks),
                disp(['Analysing subject ' num2str(sj) ', session ' num2str(session) ', block ' num2str(block)]);
                
                edffile   = dir(sprintf('P%d_s%d_b%d_*.edf', sj, session, block));
                ascfile   = dir(sprintf('P%d_s%d_b%d_*.asc', sj, session, block));
                
                if ~exist(sprintf('~/Data/RT_DDM/P%02d/S%d/P%02d_s%d_b%02d_eyeclean.mat', sj, session, sj, session, block), 'file'),
                    if ~exist(sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block), 'file'),
                        
                        % specify the filename
                        if ~exist(ascfile.name, 'file'),
                            % CONVERT TO ASC
                            if exist('~/code/Tools/eye/edf2asc-linux', 'file'),
                                system(sprintf('%s %s -input', '~/code/Tools/eye/edf2asc-linux', edffile.name));
                            else
                                system(sprintf('%s %s -input', '~/Dropbox/code/Tools/eye/edf2asc-mac', edffile.name));
                            end
                            ascfile   = dir(sprintf('P%d_s%d_b%d_*.asc', sj, session, block));
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % making a FieldTrip structure out of EyeLink data
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % read in the asc EyeLink file
                        asc = read_eyelink_ascNK_AU(ascfile.name);
                        
                        % create events and data structure
                        [data, event] = asc2dat(asc);
                        
                        % save
                        savefast(sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block), 'data', 'asc', 'event');
                    else
                        load(sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block));
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % linear blink interpolation
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    close all;
                    
                    [data, ~] = blinkinterpolate_mathot(asc, ...
                        data, 1, sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block), 1500); % padding
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % lowpass filter
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % fixdata             = data; % keep unfiltered data for fixation error rejection
                    cfg                 = [];
                    cfg.bpfilter        = 'yes';
                    cfg.bpfreq          = [0.01 4]; % highpass as 0.01 to get rid of superslow fluctuations, lowpass at 4 to get rid of fast noise
                    cfg.bpfiltord       = 2; % third order butterworth not stable
                    data                = ft_preprocessing(cfg, data);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % zscore over this block
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    pupilchan   = find(strcmp(data.label, 'EyePupil')==1);
                    pupildat    = data.trial{1}(pupilchan,:);
                    
                    % ZSCORE
                    pupiltimecourse = zscore(pupildat); % zscore
                    data.trial{1}(pupilchan,:) = pupiltimecourse; % put back in
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % define trials
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % same idea as for MEG data
                    cfg                         = [];
                    cfg.trialfun                = 'trialfun_EL_Elisa';
                    cfg.trialdef.pre            = 0;
                    cfg.trialdef.post           = 2;
                    cfg.event                   = event;
                    cfg.dataset                 = ascfile.name;
                    cfg.fsample                 = asc.fsample;
                    cfg.sj                      = sj;
                    cfg.session                 = session;
                    [cfg]                       = ft_definetrial(cfg);
                    
                    data                        = ft_redefinetrial(cfg, data); %make trials
                    data.trialinfo              = cfg.trl(:,4:end);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % mark bad trials with blinks between fix and resp
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    rejectartefacts = false;
                    if rejectartefacts,
                        
                        cfg = [];
                        cfg.artfctdef.reject                = 'nan'; % pad with nans
                        cfg.artfctdef.feedback              = 'summary'; % feedback will go into diary
                        cfg.artfctdef.blink.artifact        = blinksmp;
                        data2 = ft_rejectartifact(cfg, data);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % mark bad trials with saccades between fix and resp
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        saccthreshold = nanmean(saccsmp(:,3)) + nanstd(saccsmp(:,3)); % 69 pixels is 2 deg of visual angle (inner radius)
                        
                        cfg = [];
                        cfg.artfctdef.reject                = 'nan'; % pad with nans
                        cfg.artfctdef.feedback              = 'summary'; % feedback will go into diary
                        cfg.artfctdef.saccade.artifact      = saccsmp(find(saccsmp(:,3)>=saccthreshold), :);
                        data = ft_rejectartifact(cfg, data);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % mark bad trials that break fixation
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        fixdat              = cell2mat(fixdata.trial);
                        % often, the gaze position isn't exactly in the
                        % middle of the screen - so take the median gaze
                        % per block
                        xmean               = nanmedian(fixdat(1,:));
                        ymean               = nanmedian(fixdat(2,:));
                        xfixrange           = 3*nanstd(fixdat(1,:)); % 4 degrees of visual angle, 70 pixels (radius of inner circle around fixpoint)
                        yfixrange           = 3*nanstd(fixdat(2,:));
                        
                        for t = 1:length(data.trial),
                            clear samples2reject
                            % find samples to be padded
                            samples2reject = find(data.trial{t}(1,:) > xmean+xfixrange | ...
                                data.trial{t}(1,:) < xmean-xfixrange | ...
                                data.trial{t}(2,:) > ymean+yfixrange | ...
                                data.trial{t}(2,:) < ymean-yfixrange);
                            if ~isempty(samples2reject),
                                data.trial{t}(1:3,samples2reject) = NaN;
                                fprintf('Padding fixation break nans in trial %d \n', data.trialinfo(t,11));
                            end
                        end
                        
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % downsample before saving
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    cfg             = [];
                    cfg.resamplefs  = 50;
                    cfg.fsample     = data.fsample;
                    
                    % see Niels' message on the FT mailing list
                    samplerows = [1 3 6 8]; %find(data.trialinfo(1,:)>100); %indices of the rows with sample values (and not event codes)
                    data.trialinfo(:,samplerows) = round(data.trialinfo(:,samplerows) * (cfg.resamplefs/cfg.fsample));
                    
                    disp('resampling data...');
                    % manually downsample, to avoid interpolating spike problems
                    step = cfg.fsample/cfg.resamplefs;
                    for t = 1:length(data.trial),
                        data.trial{t} = data.trial{t}(:,1:step:end);
                        %  fixdata.trial{t} = fixdata.trial{t}(:, 1:step:end);
                        data.time{t} = data.time{t}(:,1:step:end);
                        % fixdata.time{t} = fixdata.time{t}(:,1:step:end);
                        
                    end
                    data.fsample    = cfg.resamplefs;
                    %fixdata.fsample = cfg.resamplefs;
                    
                    cd ..
                    disp(['Saving... ' sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block)]);
                    % save these datafiles before appending
                    savefast(sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block), 'data');
                    cd(['S' num2str(session)]);
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % now append all the eyelink files together
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        eyelinkfiles = dir(sprintf('~/Data/RT_DDM/P%02d/P%02d*eyeclean.mat', sj, sj));
        cd(sprintf('~/Data/RT_DDM/P%02d/', sj));
        cfg = [];
        cfg.inputfile = {eyelinkfiles.name};
        cfg.outputfile = sprintf('P%02d_alleye.mat', sj);
        ft_appenddata(cfg);
        
        % save some space
        delete *eyeclean.mat
        
        toc
    end
    cd ..
    
end

