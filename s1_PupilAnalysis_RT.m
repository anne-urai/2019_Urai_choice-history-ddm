function [] = s1_PupilAnalysis_RT(subjects)
% preprocess pupil data again
% Anne Urai, 2017

if ~isdeployed,
    addpath(genpath('~/code/RT_RDK'));
    addpath(genpath('~/code/Tools'));
    addpath('~/Documents/fieldtrip/');
    ft_defaults;
end

if ischar('subjects'), subjects = str2num(subjects); end
% path
datapath = '~/Data/RT_RDK';
if ~exist('subjects', 'var'),
    subjects =  [3     4     5     6     7     8     9    10    11    12    13   ...
        14    15    17    18    19    20    21    22    23    24    25];
end

% determine which subject this job will work on
for sj = (unique(subjects)),

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

            edffile   = dir(sprintf('%s/P%02d/pupil/P%d_s%d_b%d_*.edf', datapath, sj, sj, session, block));
            ascfile   = dir(sprintf('%s/P%02d/pupil/P%d_s%d_b%d_*.asc', datapath, sj, sj, session, block));

            % specify the filename
            if ~exist(sprintf('%s/P%02d/pupil/%s', datapath, sj, ascfile.name), 'file'),
                % CONVERT TO ASC
                if exist('~/code/Tools/eye/edf2asc-linux', 'file'),
                    system(sprintf('%s %s -input -failsafe', '~/code/Tools/eye/edf2asc-linux', edffile.name));
                else
                    system(sprintf('%s %s -input', '~/Dropbox/code/Tools/eye/edf2asc-mac', edffile.name));
                end
                ascfile   = dir(sprintf('%s/P%02d/pupil/P%d_s%d_b%d_*.asc', datapath, sj, sj, session, block));
            end

            % ==================================================================
            % making a FieldTrip structure out of EyeLink data
            % ==================================================================

            clear blinksmp saccsmp

            % read in the asc EyeLink file
            asc = read_eyelink_ascNK_AU(sprintf('%s/P%02d/pupil/%s', datapath, sj, ascfile.name));

            % create events and data structure, parse asc
            [data, event, blinksmp, saccsmp] = asc2dat(asc);

            % ==================================================================
            % blink interpolation
            % ==================================================================

            newpupil = blink_interpolate(data, blinksmp, 0);
            data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = newpupil;

            % suplabel(sprintf('P%02d-S%d-b%d', sj, session, block), 't');
            % saveas(gcf, sprintf('%s/Figures/P%02d_s%d_b%d_preproc.pdf', datapath, sj, session, block), 'pdf');

            % ==================================================================
            % regress out pupil response to blinks and saccades
            % ==================================================================

            % for this, use only EL-defined blinksamples
            % dont add back slow drift for now - baseline doesn't make a lot of
            % sense after this anymore...
            addBackSlowDrift = 0;

            pupildata = data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyepupil')),:);
            newpupil = blink_regressout(pupildata, data.fsample, blinksmp, saccsmp, 0, addBackSlowDrift);
            % put back in fieldtrip format
            data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyepupil')),:) = newpupil;
            % drawnow;
            % saveas(gcf,  sprintf('%s/Figures/P%02d_s%d_b%d_projectout.pdf', datapath, sj, session, block), 'pdf');

            % ==================================================================
            % zscore since we work with the bandpassed signal
            % ==================================================================

            data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = ...
                zscore(data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:));

            % ==================================================================
            % define trials
            % ==================================================================

            cfg                         = [];
            cfg.trialfun                = 'trialfun_EL_Elisa';
            cfg.trialdef.pre            = 0;
            cfg.trialdef.post           = 42;
            cfg.event                   = event;
            cfg.dataset                 = ascfile.name;
            cfg.fsample                 = asc.fsample;
            cfg.sj                      = sj;
            cfg.session                 = session;
            trialfun_EL_Elisa(cfg); % so that the deployed version includes this func
            [cfg]                       = ft_definetrial(cfg);

            data                        = ft_redefinetrial(cfg, data); %make trials
            data.trialinfo              = cfg.trl(:,4:end);

            % ==================================================================
            % downsample before saving
            % ==================================================================

            cfg             = [];
            cfg.resamplefs  = 100;
            cfg.fsample     = data.fsample;

            % see Niels' message on the FT mailing list
            samplerows = find(mean(data.trialinfo) > 1000); % indices of the rows with sample values (and not event codes)
            data.trialinfo(:,samplerows) = round(data.trialinfo(:,samplerows) * (cfg.resamplefs/cfg.fsample));

            % use fieldtrip to resample
            data = ft_resampledata(cfg, data);

            % save these datafiles before appending
            fprintf('Saving %s/P%02d/pupil/P%02d_s%d_b%02d_eyeclean.mat \n', datapath, sj, sj, session, block);
            savefast(sprintf('%s/P%02d/pupil/P%02d_s%d_b%02d_eyeclean.mat', datapath, sj,  sj, session, block), 'data');

        end
    end
if 0,

% ==================================================================
% now append all the eyelink files together
% ==================================================================

% check if the full dataset is not there yet
eyelinkfiles = rdir(sprintf('%s/P%02d/pupil/P%02d*_eyeclean.mat', datapath, sj, sj));

% make sure these are in the right order!
% otherwise, indexing of trials will go awry
for f = 1:length(eyelinkfiles),
    scandat         = sscanf(eyelinkfiles(f).name, 'P%*d_s%d_b%d*.mat');
    snum(f,:)       = scandat';
end
[sorted, sortidx]   = sort(snum(:,1)); % sort by session
sorted(:,2)         = snum(sortidx, 2); % sort by block
eyelinkfiles        = eyelinkfiles(sortidx);

cfg             = [];
cfg.inputfile   = {eyelinkfiles.name};
cfg.outputfile  = sprintf('%s/P%02d_alleye.mat', datapath, sj);
ft_appenddata(cfg);

% to save disk space
% eyelinkfiles = rdir(sprintf('%s/P%02d/pupil/P%02d*_eyeclean.mat', datapath, sj, sj));

end
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
end

end
