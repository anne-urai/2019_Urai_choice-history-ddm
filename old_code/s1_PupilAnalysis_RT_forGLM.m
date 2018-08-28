function [] = s1_PupilAnalysis_RT_forGLM(subjects)
% preprocess all data again, now resulting in a continuous file that I can
% use for GLM

addpath(genpath('~/Dropbox/code/RT_RDK'));
addpath('~/Documents/MATLAB/toolboxes/fieldtrip/');
ft_defaults;
clc; dbstop if error;

if ~exist('subjects', 'var'),
    % get all those subjects who have a folder
    cd ~/Data/RT_DDM/;
    s = dir('P*');
    s = {s([s(:).isdir]).name};
    for i = 1:length(s), subjects(i) = str2num(s{i}(2:3)); end
end
tic;

% determine which subject this job will work on
for sj = (unique(subjects)),
    
    cd(sprintf('~/Data/RT_DDM/P%02d/', sj));
    delete *eyeclean.mat
    
    clc;
    % ==================================================================
    % REORGANIZE DATA LOCALLY
    % ==================================================================
    sessions = 1:5;
    
    if sj == 4,
        sessions = 1:7;
    elseif  sj == 15,
        sessions = 1:5;
    elseif sj == 17,
        sessions = [1 2 3 5];
    elseif sj == 18,
        sessions = 2:5;
    elseif sj == 20,
        sessions = 2:5;
    end
    
    % ==================================================================
    % LOAD IN SUBJECT SPECIFICS AND READ DATA
    % ==================================================================
    
    for session = unique(sessions),
        
        cd(sprintf('~/Data/RT_DDM/P%02d/S%d/', sj, session));
        
        blocks = 1:10;
        
        % exceptions !!!
        if sj == 4 && session == 2,
            blocks = [1:6 9];
        elseif sj == 8 && session == 2,
            blocks = [1 2 4 5 6 7 8 9 10];
        elseif sj == 10 && session < 3,
            blocks = 2:10;
        elseif sj == 10 && session == 3,
            blocks = [1:4 6:10];
        elseif sj == 12 && session == 1,
            blocks = 1:7;
        elseif sj == 12 && session == 4,
            blocks = 1:9;
        elseif sj == 13 && session == 1,
            blocks = 1:4;
        elseif sj == 13 && session == 2,
            blocks = 2:10;
        elseif sj == 13 && session == 3,
            blocks = 2:9;
        elseif sj == 13 && session == 4,
            blocks = 1:8;
        elseif sj == 14 && session == 1,
            blocks = 1:3;
        elseif sj == 14 && session == 2,
            blocks = 1:8;
        elseif sj == 14 && session == 4,
            blocks = 2:10;
        elseif sj == 14 && session == 3,
            blocks = [1:5 7:10];
        elseif sj == 15 && session == 1,
            blocks = 1:5;
        elseif sj == 15 && session == 5,
            blocks = 1:8;
        elseif sj == 17 && session == 1,
            blocks = 1:7;
        elseif sj == 19 && session == 1,
            blocks = 1:4;
        elseif sj == 21 && session == 1,
            blocks = 1:6;
        elseif sj == 22 && session == 1,
            blocks = 1:6;
        elseif sj == 23 && session == 1,
            blocks = 1:4;
        elseif sj == 24 && session == 1,
            blocks = 1:5;
        elseif sj == 25 && session == 1,
            blocks = 1:7;
        end
        
        for block = unique(blocks),
            disp(['Analysing subject ' num2str(sj) ', session ' num2str(session) ', block ' num2str(block)]);
            
            edffile   = dir(sprintf('P%d_s%d_b%d_*.edf', sj, session, block));
            ascfile   = dir(sprintf('P%d_s%d_b%d_*.asc', sj, session, block));
            
            if ~exist(sprintf('~/Data/RT_DDM/P%02d/P%02d_s%d_b%02d_eyeclean.mat', sj, sj, session, block), 'file'),
                if ~exist(sprintf('~/Data/RT_DDM/P%02d/S%d/P%02d_s%d_b%02d_eye.mat', sj, session, sj, session, block), 'file'),
                    
                    % specify the filename
                    if ~exist(ascfile.name, 'file'),
                        % CONVERT TO ASC
                        if exist('~/code/Tools/eye/edf2asc-linux', 'file'),
                            % on UKE cluster
                            system(sprintf('%s %s -input', '~/code/Tools/eye/edf2asc-linux', edffile.name));
                        else
                            % on local MacbookPro
                            system(sprintf('%s %s -input', '~/Dropbox/code/Tools/eye/edf2asc-mac', edffile.name));
                        end
                        ascfile   = dir(sprintf('~/Data/RT_DDM/P%02d/S%d/P%d_s%d_b%d_*.asc', sj, session, sj, session, block));
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % making a FieldTrip structure out of EyeLink data
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % read in the asc EyeLink file
                    asc = read_eyelink_ascNK_AU(ascfile.name);
                    
                    % create events and data structure
                    [data, event] = asc2dat(asc);
                    
                    % save
                    savefast(sprintf('~/Data/RT_DDM/P%02d/S%d/P%02d_s%d_b%02d_eye.mat', sj, session, sj, session, block), 'data', 'asc', 'event');
                else
                    load(sprintf('~/Data/RT_DDM/P%02d/S%d/P%02d_s%d_b%02d_eye.mat', sj, session, sj, session, block));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % linear blink interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                close all;
                [data, ~] = blinkinterpolate_mathot(asc, ...
                    data, 1, sprintf('P%02d_s%d_b%02d_eye.mat', sj, session, block), 1000);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % bandpass filter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % only keep the pupil
                cfg                 = [];
                cfg.channel         = 'EyePupil';
                data                = ft_selectdata(cfg, data);
                
                cfg                 = [];
                cfg.bpfilter        = 'yes';
                cfg.bpfreq          = [0.05 4]; % as in JWs paper
                cfg.bpfiltord       = 2;
                data                = ft_preprocessing(cfg, data);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % zscore
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                data.trial{1} = zscore(data.trial{1});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % define trials
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % same idea as for MEG data
                cfg                         = [];
                cfg.trialfun                = 'trialfun_EL_Elisa';
                cfg.trialdef.pre            = 0; % start at fixation
                cfg.trialdef.post           = 2; % 2s after feedback
                cfg.event                   = event;
                cfg.dataset                 = ascfile.name;
                cfg.fsample                 = asc.fsample;
                cfg.sj                      = sj;
                cfg.session                 = session;
                [cfg]                       = ft_definetrial(cfg);
                trialinfo                   = cfg.trl; % save for later
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % downsample before saving
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                cfg             = [];
                cfg.resamplefs  = 100; % the larger the sampling rate, the slower deconvolution gets
                cfg.fsample     = data.fsample;
                cfg.detrend     = 'no';
                cfg.demean      = 'yes'; % will zscore per block anyway
                
                % see Niels' message on the FT mailing list
                samplerows              = [1 3 7 9]; % these have samples in them
                trialinfo(:,samplerows) = round(trialinfo(:,samplerows) * (cfg.resamplefs/cfg.fsample));
                blinksmp                = round(blinksmp .* (cfg.resamplefs/cfg.fsample));
                
                disp('resampling data...');
                data           = ft_resampledata(cfg, data);
                data.trialinfo = trialinfo;
                data.blinksmp  = blinksmp;
                
                disp(['Saving... ' sprintf('P%02d-S%02d_b%02d_eyeclean_continuous.mat', sj, session, block)]);
                % save these datafiles before appending
                savefast(sprintf('~/Data/RT_DDM/P%02d/P%02d-S%02d_b%02d_eyeclean_continuous.mat', sj, sj, session, block), 'data');
                
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now append all the eyelink files together
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clc;
    clearvars -except sj subjects
    
    cd(sprintf('~/Data/RT_DDM/P%02d', sj));
    eyelinkfiles        = dir('*_eyeclean_continuous.mat');
    % make sure these are in the right order!
    % otherwise, indexing of trials will go awry
    for f = 1:length(eyelinkfiles),
        scandat         = sscanf(eyelinkfiles(f).name, 'P%*d-S%d_b%d*.mat');
        snum(f,:)       = scandat';
    end
    [sorted, sortidx]   = sort(snum(:,1)); % sort by session
    sorted(:,2)         = snum(sortidx, 2); % sort by block
    eyelinkfiles        = eyelinkfiles(sortidx);
    
    cfg = [];
    cfg.inputfile   = {eyelinkfiles.name};
    alldata         = ft_appenddata(cfg);
    
    % also get a concatenated trialinfo
    alltrl = []; allblink = [];
    toAdd  = 0;
    samplerows = [1 3 7 9]; % these have samples in them
    
    for f = 1:length(eyelinkfiles),
        clear data;
        load(eyelinkfiles(f).name);
        trl         = data.trialinfo;
        blinksmp    = data.blinksmp;
        
        % add the length of the session before
        trl(:, samplerows)  = trl(:, samplerows) + toAdd;
        blinksmp            = blinksmp + toAdd;
        
        toAdd               = toAdd + numel(data.time{1});
        alltrl              = [alltrl; trl];
        allblink            = [allblink; blinksmp];
    end
    
    % remove weird stuff
    toremove = find(isnan(alltrl(:, 1)));
    alltrl(toremove, :) = [];
    
    % reappend
    alldata.trialinfo          = alltrl;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % append trialinfo that can be used to epoch later
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data           = alldata;
    data.trialinfo = alltrl;
    data.blinksmp  = allblink;
    
    trial   = cat(2, data.trial{:});
    time    = cat(2, data.time{:});
    data    = rmfield(data, {'trial', 'time'});
    data.trial{1}   = trial;
    data.time{1}    = time;
    
    savefast(sprintf('~/Data/RT_DDM/P%02d_alleye_continuous.mat', sj), 'data');
    fprintf('P%02d_alleye_continuous.mat \n', sj);
    
    toc;
end
end
