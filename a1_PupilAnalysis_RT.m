function [] = s1_PupilAnalysis_RT(sj, session, block)
  % preprocess pupil data again
  % Anne Urai, 2017

  if ~isdeployed,
    addpath(genpath('~/code/RT_RDK'));
    addpath(genpath('~/code/Tools'));
    addpath('~/Documents/fieldtrip/');
    ft_defaults;
  end

  tic;
  if ischar(sj),      sj = str2num(sj);
  end
  if ischar(session), session = str2num(session);
  end
  if ischar(block),   block = str2num(block);
  end

  % path
  datapath = '~/Data/RT_RDK';

  % ==================================================================
  % LOAD IN SUBJECT SPECIFICS AND READ DATA
  % ==================================================================

  disp(['Analysing subject ' num2str(sj) ', session ' num2str(session) ', block ' num2str(block)]);

  edffile   = dir(sprintf('%s/P%02d/pupil/P%d_s%d_b%d_*.edf', datapath, sj, sj, session, block));
  ascfile   = dir(sprintf('%s/P%02d/pupil/P%d_s%d_b%d_*.asc', datapath, sj, sj, session, block));

  % if length(edffile) > 1, edffile = edffile(1); end % fix later
  % don't run this script if these files dont exist
  if (isempty(edffile) && isempty(ascfile)),
    disp('skipping');
    return;
  end

  % specify the filename
  if ~isempty(edffile),
    % CONVERT TO ASC
    disp('converting EDF2ASC');
    system(sprintf('%s %s -input -failsafe -y', '~/code/Tools/eye/edf2asc-linux', ...
    sprintf('%s/P%02d/pupil/%s', datapath, sj, edffile.name)));
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

  %suplabel(sprintf('P%02d-S%d-b%d', sj, session, block), 't');
  %saveas(gcf, sprintf('%s/Figures/P%02d_s%d_b%d_preproc.pdf', datapath, sj, session, block), 'pdf');

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
  %drawnow;
  %saveas(gcf,  sprintf('%s/Figures/P%02d_s%d_b%d_projectout.pdf', datapath, sj, session, block), 'pdf');

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
  toc;
end
