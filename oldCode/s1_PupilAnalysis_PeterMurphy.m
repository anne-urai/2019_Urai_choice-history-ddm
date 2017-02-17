function [] = s1_PupilAnalysis_PeterMurphy()
% instead of normalizing per block, do this over one whole session to
% approximate the pipeline that Peter Murphy used

addpath(genpath('~/Dropbox/code/RT_RDK'));
addpath('~/Documents/MATLAB/toolboxes/fieldtrip/');
ft_defaults;
clc; dbstop if error;

if ~exist('subjects', 'var'),
    % get all those subjects who have a folder
    cd ~/Data/PeterMurphy/;
    s = dir('P*');
    s = {s([s(:).isdir]).name};
    for i = 1:length(s), subjects(i) = str2num(s{i}(2:4)); end
end


for sj = unique(subjects),
    
    cd(sprintf('~/Data/PeterMurphy/P%02d/', sj));
    
    % ==================================================================
    % LOAD IN SUBJECT SPECIFICS AND READ DATA
    % ==================================================================
    
    blocks = 1:5;
    
    for block = unique(blocks),
        disp(['Analysing subject ' num2str(sj) ', block ' num2str(block)]);
        
        edffile   = dir(sprintf('%d_%d.edf', sj, block));
        ascfile   = dir(sprintf('%d_%d.asc', sj, block));
        
        if ~exist(sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eyeclean.mat', sj, sj, block), 'file'),
            if ~exist(sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eye.mat', sj, sj, block), 'file'),
                
                % specify the filename
                if ~exist(ascfile.name, 'file'),
                    % CONVERT TO ASC
                    if exist('~/code/Tools/eye/edf2asc-linux', 'file'),
                        system(sprintf('%s %s -input', '~/code/Tools/eye/edf2asc-linux', edffile.name));
                    else
                        system(sprintf('%s %s -input', '~/Dropbox/code/Tools/eye/edf2asc-mac', edffile.name));
                    end
                    ascfile   = dir(sprintf('~/Data/PeterMurphy/P%02d/%d_%d.asc', sj, sj, block));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % making a FieldTrip structure out of EyeLink data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % read in the asc EyeLink file
                asc = read_eyelink_ascNK_AU(ascfile.name);
                
                % create events and data structure
                [data, event] = asc2dat(asc);
                
                % save
                savefast(sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eye.mat', sj, sj,  block), 'data', 'asc', 'event');
            else
                load(sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eye.mat', sj, sj, block));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % linear blink interpolation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            close all;
            [data, ~] = blinkinterpolate_mathot(asc, ...
                data, 1, sprintf('P%02d_b%02d_eye.mat', sj, block), 1000);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % lowpass filter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            cfg                 = [];
            cfg.lpfilter        = 'yes';
            cfg.lpfreq          = 10; % highpass as 0.01 to get rid of superslow fluctuations, lowpass at 4 to get rid of fast noise
            datafilt            = ft_preprocessing(cfg, data);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % define trials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % same idea as for MEG data
            cfg                         = [];
            cfg.trialfun                = 'trialfun_EL_PeterMurphy';
            cfg.trialdef.pre            = 1; % from 1s before stim onset
            cfg.trialdef.post           = 2; % to 2s after response
            cfg.event                   = event;
            cfg.dataset                 = ascfile.name;
            cfg.fsample                 = asc.fsample;
            cfg.sj                      = sj;
            cfg.block                   = block;
            [cfg]                       = ft_definetrial(cfg);
            
            data                        = ft_redefinetrial(cfg, datafilt); %make trials
            data.trialinfo              = cfg.trl(:,4:end);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % downsample before saving
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if 0,
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
            end
            
            disp(['Saving... ' sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eyeclean.mat', sj, sj, block)]);
            % save these datafiles before appending
            savefast(sprintf('~/Data/PeterMurphy/P%02d/P%02d_b%02d_eyeclean.mat', sj, sj, block), 'data');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % at the end of each session, look at the full timecourse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd(sprintf('~/Data/PeterMurphy/P%02d/', sj));
    eyelinkfiles = dir(sprintf('P%02d_b*eyeclean.mat', sj));
    
    % make sure these are in the right order!
    % otherwise, indexing of trials will go awry
    clear snum
    for f = 1:length(eyelinkfiles),
        scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_b%d*.mat');
        snum(f,:)       = scandat';
    end
    [sorted, sortidx]   = sort(snum(:,1)); % sort by session
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
    % compute percent signal change over all data from this SJ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('transforming into % signal change ...');
    medianpupil = median(cat(2, data.trial{:}));
    
    for t = 1:length(data.trial),
        data.trial{t} = (data.trial{t} - medianpupil) ./ medianpupil * 100; % put back in
    end
    
    % visualize the signal
    subplot(512);
    plot(cat(2, data.trial{:}), '.'); axis tight; box off;
    
    disp(['Saving... ' sprintf('P%02d_alleye.mat', sj)]);
    % save these datafiles before appending
    savefast(sprintf('P%02d_alleye.mat', sj), 'data');
    
    % move upwards
    movefile(sprintf('~/Data/PeterMurphy/P%02d/P%02d_alleye.mat', sj, sj), sprintf('~/Data/PeterMurphy/P%02d_alleye.mat', sj));
    
end


end
