function [] = s1_PupilAnalysis_RT(subjects)
% instead of normalizing per block, do this over one whole session to
% approximate the pipeline that Peter Murphy used

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
subjects(subjects<23) = [];
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
                            system(sprintf('%s %s -input', '~/code/Tools/eye/edf2asc-linux', edffile.name));
                        else
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
                % lowpass filter
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % fixdata             = data; % keep unfiltered data for fixation error rejection
                cfg                 = [];
                cfg.lpfilter        = 'yes';
                cfg.lpfreq          = 10; % highpass as 0.01 to get rid of superslow fluctuations, lowpass at 4 to get rid of fast noise
                datafilt            = ft_preprocessing(cfg, data);
                
                clf; subplot(511);
                pupilchan   = find(strcmp(datafilt.label, 'EyePupil')==1);
                plot(data.time{1}, data.trial{1}(pupilchan, :));
                hold on; plot(datafilt.time{1}, datafilt.trial{1}(pupilchan, :));
                axis tight; box off; drawnow;
                
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
                
                data                        = ft_redefinetrial(cfg, datafilt); %make trials
                data.trialinfo              = cfg.trl(:,4:end);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % downsample before saving
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                cfg             = [];
                cfg.resamplefs  = 100;
                cfg.fsample     = data.fsample;
                
                % see Niels' message on the FT mailing list
                samplerows = find(data.trialinfo(1,:)>100); %indices of the rows with sample values (and not event codes)
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
                
                cd ..
                disp(['Saving... ' sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block)]);
                % save these datafiles before appending
                savefast(sprintf('P%02d_s%d_b%02d_eyeclean.mat', sj, session, block), 'data');
                cd(['S' num2str(session)]);
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % at the end of each session, look at the full timecourse
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cd(sprintf('~/Data/RT_DDM/P%02d/', sj));
        eyelinkfiles = dir(sprintf('P%02d_s%d*_b*eyeclean.mat', sj, session));
        
        % make sure these are in the right order!
        % otherwise, indexing of trials will go awry
        clear snum
        for f = 1:length(eyelinkfiles),
            scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_s%d_b%d*.mat');
            snum(f,:)       = scandat';
        end
        [sorted, sortidx]   = sort(snum(:,1)); % sort by session
        sorted(:,2)         = snum(sortidx, 2); % sort by block
        eyelinkfiles        = eyelinkfiles(sortidx);
        
        cfg                 = [];
        cfg.inputfile       = {eyelinkfiles.name};
        data                = ft_appenddata(cfg);
        pupilchan           = find(strcmp(data.label, 'EyePupil')==1);
        
        % keep only the pupil
        cfg                 = [];
        cfg.channel         = 'EyePupil';
        data                = ft_selectdata(cfg, data);
        
        clf;
        subplot(511);
        plot(cat(2, data.trial{:}), '.'); axis tight; box off;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute percent signal change over the whole session
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('transforming into % signal change ...');
        medianpupil = median(cat(2, data.trial{:}));
        
        for t = 1:length(data.trial),
            data.trial{t} = (data.trial{t} - medianpupil) ./ medianpupil * 100; % put back in
        end
        
        % visualize the signal
        subplot(512);
        plot(cat(2, data.trial{:}), '.'); axis tight; box off;
        % saveas(gcf, sprintf('~/Data/RT_DDM/P%02d/P%02d_s%d_allpreproc.png', sj, sj, session), 'png');
        
        disp(['Saving... ' sprintf('P%02d_s%d_allblocks_eyeclean.mat', sj, session)]);
        % save these datafiles before appending
        savefast(sprintf('P%02d_s%d_allblocks_eyeclean.mat', sj, session), 'data');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now append all the eyelink files together
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % check if the full dataset is not there yet
    clear eyelinkfiles snum
    
    cd(sprintf('~/Data/RT_DDM/P%02d/', sj));
    eyelinkfiles = dir(sprintf('P%02d*_allblocks_eyeclean.mat', sj));
    
    % make sure these are in the right order!
    % otherwise, indexing of trials will go awry
    for f = 1:length(eyelinkfiles),
        scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_s%d*.mat');
        snum(f,:)       = scandat';
    end
    [sorted, sortidx]   = sort(snum(:,1)); % sort by session
    eyelinkfiles        = eyelinkfiles(sortidx);
    
    cfg = [];
    cfg.inputfile = {eyelinkfiles.name};
    cfg.outputfile = sprintf('P%02d_alleye.mat', sj);
    ft_appenddata(cfg);

    % move upwards
    movefile(sprintf('~/Data/RT_DDM/P%02d/P%02d_alleye.mat', sj, sj), sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CORRECT SESSION NRS FOR SOME PARTICIPANTS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clearvars -except sj subjects
    load(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj));
    
    sessionsdone = unique(data.trialinfo(:, end));
    switch sj
        case 4
            assert(isequal(sessionsdone', 1:7), 'did not find 5 sessions');
        case 17
            data.trialinfo(find(data.trialinfo(:, end)==5), end) = 4;
            sessionsdone = unique(data.trialinfo(:, end));
            assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
        case 18
            data.trialinfo(:, end) = data.trialinfo(:, end) - 1;
            sessionsdone = unique(data.trialinfo(:, end));
            assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
        case 20
            data.trialinfo(:, end) = data.trialinfo(:, end) - 1;
            sessionsdone = unique(data.trialinfo(:, end));
            assert(isequal(sessionsdone', 1:4), 'did not find 5 sessions');
        otherwise
            assert(isequal(sessionsdone', 1:5), 'did not find 5 sessions');
    end
    
    save(sprintf('~/Data/RT_DDM/P%02d_alleye.mat', sj), 'data');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT ALL
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data.fsample    = 100;
    pupilchan       = find(strcmp(data.label, 'EyePupil')==1);
    
    clf;
    for sess = unique(data.trialinfo(:, end))',
        dat = cat(2, data.trial{data.trialinfo(:,end)==sess});
        subplot(max(unique(data.trialinfo(:, end))), 1, sess); plot(dat(pupilchan, :), '.');
        axis tight; box off; set(gca, 'tickdir', 'out');
        drawnow;
    end
    suplabel(sprintf('P%02d', sj), 't');
    suplabel('Percent signal change', 'y');
    
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, '-dpdf', '-painters', sprintf('~/Data/RT_DDM/Figures/P%02d_AllPupilCheck.pdf', sj));
    
end
end
