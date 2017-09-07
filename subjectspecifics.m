function subjectdata = subjectspecifics(sj)

  % path = '/projects/0/neurodec/Data/MEG-PL';
  % for sj = 2:65,
  %   files = dir(sprintf('/projects/0/neurodec/home/aurai/Data/MEG-PL/P%02d/MEG/Preproc/*_data.mat', sj));
  %   for f = 1:length(files),
  %     movefile(sprintf('/projects/0/neurodec/home/aurai/Data/MEG-PL/P%02d/MEG/Preproc/%s', sj, files(f).name), ...
  %     sprintf('/projects/0/neurodec/Data/MEG-PL/P%02d/MEG/Preproc/%s', sj, files(f).name));
  %   end
  % end

  % determine the path to the data
  usr = getenv('USER');
  switch usr
  case 'aurai' % uke cluster
    path = '~/Data/MEG-PL';
  case 'aeurai' % cartesius/lisa
    path = '/projects/0/neurodec/Data/MEG-PL';
  case 'anne' % macbook pro
    path = '/Users/anne/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/neurodec/Data/MEG-PL';
    path = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL'; % copied to local!
  otherwise
    error('could not find the data path');
  end
  subjectdata.path = path;

  % folder structure is identical across subjects
  subjectdata.rawdir              = sprintf('%s/P%02d/MEG/Raw', path, sj);
  subjectdata.preprocdir          = sprintf('%s/P%02d/MEG/Preproc', path, sj);
  subjectdata.lockdir             = sprintf('%s/P%02d/MEG/Locked',  path, sj);
  subjectdata.tfrdir              = sprintf('%s/P%02d/MEG/TFR',  path, sj);
  subjectdata.mridir              = sprintf('%s/P%02d/MRI',  path, sj);
  subjectdata.sourcedir           = sprintf('%s/P%02d/MEG/Source',  path, sj);
  subjectdata.eyedir              = sprintf('%s/P%02d/Eye',  path, sj);
  subjectdata.behavdir            = sprintf('%s/P%02d/Behav',  path, sj);
  subjectdata.rsdir               = sprintf('%s/P%02d/MEG/RestingState',  path, sj);
  subjectdata.dotsdir             = sprintf('%s/P%02d/DotCoord',  path, sj);
  subjectdata.statsdir            = sprintf('%s/P%02d/Stats',  path, sj);
  subjectdata.roidir              = sprintf('%s/P%02d/MEG/ROI',  path, sj);
  subjectdata.ldadir              = sprintf('%s/P%02d/MEG/LDA',  path, sj);

  % general path
  subjectdata.figsdir             = sprintf('%s/PreprocFigs', path);

  % =========================================== %
  % GRAND AVERAGES vs INDIVIDUAL
  % =========================================== %

  if ~isnumeric(sj),

    % drug groups
    subjectdata.all         = [2:65]; % exclude P16, 11 and 37
    subjectdata.all(ismember(subjectdata.all, [11 37 16])) = [];
    subjectdata.placebo     = [3,6,9,12,15,18,21,23,26,29,32,35,38,41,44,47,50,56,59];
    subjectdata.atomoxetine = [2,5,8,14,17,20,24,27,30,33,36,39,42,45,48,51,54,57,60,61,63,65];
    subjectdata.donepezil   = [4,7,10,13,19,22,25,28,31,34,40,43,46,49,52,53,55,58,62,64]; % exclude P16

    % to check hand counterbalancing
    subjectdata.lefthand    = intersect(2:2:65, subjectdata.all);
    subjectdata.righthand   = intersect(3:2:65, subjectdata.all);

    % serial choice bias groups - based on group seletion in E0b,
    % p(repeat) across all 5 sessions
    subjectdata.alternators = [2,3,5,7,8,9,10,12,13,14,18,20,21,22,24,28,29,30,31,33,34,40,42,43,46,50,58,61];
    subjectdata.repeaters   = [4,6,15,17,19,23,25,26,27,32,35,36,38,39,41,44,45,47,48,49,51,52,53,54,55,56,57,59,60,62,63,64,65];
    assert(isequal(sort(subjectdata.all), sort([subjectdata.alternators subjectdata.repeaters])), 'subject nr mismatch');
    
    % serial choice bias only defined on MEG sessions
    subjectdata.repeatersM = [2,3,4,5,6,7,10,13,14,15,17,19,21,22,23,24,25,26,27,28,32,34,35,36,38,41,42,44,45,46,47,48,49,51,52,53,54,55,56,57,59,60,61,62,63,64,65];
    subjectdata.alternatorsM = [8,9,12,18,20,29,30,33,39,40,43,50,58];

    subjectdata.tfrdir      = sprintf('%s/GrandAverage/TFR', path);
    subjectdata.lockdir     = sprintf('%s/GrandAverage/Locked', path);
    subjectdata.statsdir    = sprintf('%s/GrandAverage/Stats',  path);
    subjectdata.sourcedir   = sprintf('%s/GrandAverage/Source',  path);
    subjectdata.mridir      = sprintf('%s/GrandAverage/MRI',  path);
    subjectdata.toidir      = sprintf('%s/GrandAverage/TOI',  path);
    subjectdata.roidir      = sprintf('%s/GrandAverage/ROI',  path);
    subjectdata.figsdir     = sprintf('%s/Figures', path);

  else
    switch lower(sj)

    case 1
      % excluded

    case 2

      subjectdata.name                = 'LB';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3]; % order for task recordings
      subjectdata.session(1).rec(2).blocks        = 1:2; % no EL on first day
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:6;
      subjectdata.session(2).rec(3).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 3

      subjectdata.name                = 'PB';
      subjectdata.drug                = 'placebo';

      subjectdata.session(1).recsorder            = 1;
      subjectdata.session(1).rec(1).blocks        = 1:10; % no EL on first day
      subjectdata.session(1).rec(1).blockidx      = [1 2 1:8]; % these were entered under a different idx
      subjectdata.session(1).rec(1).block(10).missingtrials = 51:60;
      subjectdata.session(1).restingstate         = [];

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 4

      subjectdata.name                = 'MS';
      subjectdata.drug                = 'donepezil';

      subjectdata.session(1).recsorder            = [1 2 3];
      subjectdata.session(1).rec(1).blocks        = 1;
      subjectdata.session(1).rec(2).blocks        = [2 3 4];
      subjectdata.session(1).rec(2).block(2).missingtrials = 1:12;
      subjectdata.session(1).rec(3).blocks        = [5:10];
      subjectdata.session(1).restingstate         = [];

      subjectdata.session(2).recsorder            = [1 2];
      subjectdata.session(2).rec(1).blocks        = 1:5;
      subjectdata.session(2).rec(2).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 3;

    case 5

      subjectdata.name                = 'SR';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = [1 2];
      subjectdata.session(1).rec(1).blocks        = 1:5; % no EL on first day
      subjectdata.session(1).rec(2).blocks        = 6:10;
      subjectdata.session(1).restingstate         = [];

      subjectdata.session(2).recsorder            = [1 2];
      subjectdata.session(2).rec(1).blocks        = 1:5;
      subjectdata.session(2).rec(2).blocks        = 6:10;
      subjectdata.session(2).rec(2).block(10).missingtrials = 52:60;
      subjectdata.session(2).restingstate         = 3;

    case 6

      subjectdata.name                = 'CK';
      subjectdata.drug                = 'placebo';

      subjectdata.session(1).recsorder            = 2;
      subjectdata.session(1).rec(2).blocks        = 1:5; % no EL on first day
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 7

      subjectdata.name                = 'SG';
      subjectdata.drug                = 'donepezil';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).rec(3).block(6).missingtrials = 1;
      subjectdata.session(2).restingstate         = 1;

    case 8

      subjectdata.name                = 'LL';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3];
      subjectdata.session(1).rec(2).blocks        = [1 21 22 3 4 5];
      subjectdata.session(1).rec(2).block(21).missingtrials = 25:60;
      subjectdata.session(1).rec(2).block(22).missingtrials = 42:60;

      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).rec(3).block(6).missingtrials = 1:17; % first 8 trials of block 6 missing';
      subjectdata.session(1).missingtrials        = 'block 1-5, no audio';
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 9

      subjectdata.name                = 'MM';
      subjectdata.drug                = 'placebo';

      subjectdata.session(1).recsorder            = [2 3];
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 10

      subjectdata.name                = 'HS';
      subjectdata.drug                = 'donepezil';

      subjectdata.session(1).recsorder            = [2 3];
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 11

      subjectdata.name                = 'KM';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).rec(2).block(2).missingtrials = 33;
      subjectdata.session(1).missingtrials        = 'audio off for first 3 trials of block 1';
      subjectdata.session(1).restingstate         = 1;

      % did not complete the experiment, no session 2

    case 12

      subjectdata.name                = 'LP';
      subjectdata.drug                = 'placebo';

      subjectdata.session(1).recsorder            = [4 5];
      subjectdata.session(1).rec(4).blocks        = 1:4; % rec2 and 3 dontuse, performance too low
      subjectdata.session(1).rec(5).blocks        = 5:8;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:2;
      subjectdata.session(2).rec(3).blocks        = 3:6;
      subjectdata.session(2).rec(4).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 13

      subjectdata.name                = 'MR';
      subjectdata.drug                = 'donepezil';

      subjectdata.session(1).recsorder            = 2:6;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:4;
      subjectdata.session(1).rec(4).blocks        = 5:6;
      subjectdata.session(1).rec(5).blocks        = 7:8;
      subjectdata.session(1).rec(6).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 7:9;
      subjectdata.session(2).rec(3).blocks        = [10, 1:2];
      subjectdata.session(2).rec(4).blocks        = 3:6;
      subjectdata.session(2).restingstate         = 1;

    case 14

      subjectdata.name                 = 'HH';
      subjectdata.drug                 = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3];
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).missingtrials        = 'block 1, no audio';
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:6;
      subjectdata.session(2).rec(3).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 15

      subjectdata.name                 = 'SB';
      subjectdata.drug                 = 'placebo';

      subjectdata.session(1).recsorder            = [3 4 5];
      subjectdata.session(1).rec(3).blocks        = 1:5;
      subjectdata.session(1).rec(4).blocks        = 6:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:7;
      subjectdata.session(2).rec(3).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 16

      % weird subject with loads of movement and artefacts
      % throw out!
      subjectdata.name                 = 'BB';
      subjectdata.drug                 = 'donepezil';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:6;
      subjectdata.session(1).rec(4).blocks        = 7:10;
      subjectdata.session(1).missingtrials        = 'block 3, buttons wrong way around';
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 17

      subjectdata.name                 = 'MG';
      subjectdata.drug                 = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).rec(4).blocks        = 6:7;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 18

      subjectdata.name                 = 'WF';
      subjectdata.drug                 = 'placebo';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 19

      subjectdata.name                 = 'NS';
      subjectdata.drug                 = 'donepezil';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 20

      subjectdata.name                 = 'HR';
      subjectdata.drug                 = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 21

      subjectdata.name                 = 'AE';
      subjectdata.drug                 = 'placebo';

      subjectdata.session(1).recsorder            = 6;
      subjectdata.session(1).rec(6).blocks        = 7:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 22

      subjectdata.name                  = 'MT';
      subjectdata.drug                  = 'donepezil';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 23

      subjectdata.name                  = 'FG';
      subjectdata.drug                  = 'placebo';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 24

      subjectdata.name                = 'SP';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6;
      subjectdata.session(1).rec(4).blocks        = 7:10;
      subjectdata.session(1).missingtrials        = 'edf file block 6 not complete';
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [1 2 3];
      subjectdata.session(2).rec(1).blocks        = 1:5;
      subjectdata.session(2).rec(2).blocks        = 6:9;
      subjectdata.session(2).rec(3).blocks        = 10;
      subjectdata.session(2).restingstate         = 4;

    case 25

      subjectdata.name                = 'LD';
      subjectdata.drug                = 'donepezil';

      subjectdata.session(1).recsorder            = [2 3 4];
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3];
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 26

      subjectdata.name                = 'CK';
      subjectdata.drug                = 'placebo';

      subjectdata.session(1).recsorder            = [2 3 4 5];
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).rec(4).blocks        = 6:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4];
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 27

      subjectdata.name                = 'TA';
      subjectdata.drug                = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:7;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3;
      subjectdata.session(1).rec(4).blocks        = 4:5;
      subjectdata.session(1).rec(5).blocks        = 6;
      subjectdata.session(1).rec(6).blocks        = 7:8;
      subjectdata.session(1).rec(7).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:9;
      subjectdata.session(2).rec(5).blocks        = 10;
      subjectdata.session(2).restingstate         = 1;

    case 28

      subjectdata.name                    = 'MB';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:6;
      subjectdata.session(1).rec(2).blocks        = 1;
      subjectdata.session(1).rec(3).blocks        = 2;
      subjectdata.session(1).rec(4).blocks        = 3:5;
      subjectdata.session(1).rec(5).blocks        = 6:8;
      subjectdata.session(1).rec(6).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4;
      subjectdata.session(2).rec(4).blocks        = 5:7;
      subjectdata.session(2).rec(5).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 29

      subjectdata.name                    = 'TS';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 30

      subjectdata.name                    = 'SW';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).rec(4).blocks        = 6:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:6;
      subjectdata.session(2).rec(2).blocks        = 7:10;
      subjectdata.session(2).rec(3).blocks        = 1:3;
      subjectdata.session(2).rec(4).blocks        = 4;
      subjectdata.session(2).rec(5).blocks        = 5;
      subjectdata.session(2).rec(6).blocks        = 6;
      subjectdata.session(2).restingstate         = 1;

    case 31

      subjectdata.name                    = 'KR';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 32

      subjectdata.name                    = 'LZ';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 3:6;
      subjectdata.session(1).rec(3).blocks        = 1:2;
      subjectdata.session(1).rec(4).blocks        = 3:5;
      subjectdata.session(1).rec(5).blocks        = 6:8;
      subjectdata.session(1).rec(6).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 33

      subjectdata.name 					= 'LB';
      subjectdata.drug 					= 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).rec(4).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 34

      subjectdata.name                    = 'OW';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 35

      subjectdata.name                    = 'MA';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = [6 7];
      subjectdata.session(2).rec(4).blocks        = [81 82 9 10];
      subjectdata.session(2).rec(4).block(81).missingtrials = 45:60;
      subjectdata.session(2).rec(4).block(82).missingtrials = 1:44;
      subjectdata.session(2).missingtrials        = 'in block 8, script crashed at trial 45 (2 resp buttons at the same time). Restarted once but made an error, then restarted again.';
      subjectdata.session(2).restingstate         = 1;

    case 36

      subjectdata.name                    = 'MkA';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:6;
      subjectdata.session(1).rec(4).blocks        = 7:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 37

      subjectdata.name                    = 'AE';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      % did not do session 2

    case 38

      subjectdata.name                    = 'YP';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:2;
      subjectdata.session(2).rec(3).blocks        = 3:5;
      subjectdata.session(2).rec(4).blocks        = 6:10;
      %subjectdata.session(2).rec(5).blocks        = 10;
      subjectdata.session(2).restingstate         = 1;

    case 39

      subjectdata.name                    = 'NF';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:5;
      subjectdata.session(1).rec(4).blocks        = 6:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 40

      subjectdata.name                    = 'LK';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 41

      subjectdata.name                    = 'AL';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 42

      subjectdata.name                    = 'FA';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 43

      subjectdata.name                    = 'DR';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:5;
      subjectdata.session(1).rec(4).blocks        = 6:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 44

      subjectdata.name                    = 'MS';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:4;
      subjectdata.session(1).rec(4).blocks        = 5:7;
      subjectdata.session(1).rec(5).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 45

      subjectdata.name                    = 'JR';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 46

      subjectdata.name                    = 'AH';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 47

      subjectdata.name                    = 'JR';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:6;
      subjectdata.session(2).rec(3).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 48

      subjectdata.name                    = 'JH';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 49

      subjectdata.name                    = 'AM';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:10;
      subjectdata.session(2).restingstate         = 1;

    case 50

      subjectdata.name                    = 'SS';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 51

      subjectdata.name                    = 'MS';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 52

      subjectdata.name                    = 'MH';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1;
      subjectdata.session(1).rec(3).blocks        = 2:5;
      subjectdata.session(1).rec(4).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:7;
      subjectdata.session(2).rec(3).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 53

      subjectdata.name                    = 'ACW';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 54

      subjectdata.name                    = 'MS';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = [2 3 4 5 6];
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:5;
      subjectdata.session(1).rec(4).blocks        = 6:7;
      subjectdata.session(1).rec(4).block(7).missingtrials = 58:60;
      subjectdata.session(1).rec(5).blocks        = 81;
      subjectdata.session(1).rec(5).block(81).missingtrials = 19:60;
      subjectdata.session(1).rec(6).blocks        = [82 9 10];
      %subjectdata.session(1).rec(5).block(82).missingtrials = 18:60;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:10;
      subjectdata.session(2).restingstate         = 1;

    case 55

      subjectdata.name                    = 'EI';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 3:5;
      subjectdata.session(2).rec(3).blocks        = 1:5;
      subjectdata.session(2).rec(4).blocks        = 6:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1; % had her bra on, horrible signal!

    case 56

      subjectdata.name                    = 'JH';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:4;
      subjectdata.session(1).rec(4).blocks        = 5:7;
      subjectdata.session(1).rec(5).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:8;
      subjectdata.session(2).rec(4).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 57

      subjectdata.name                    = 'HE';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:5;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:6;
      subjectdata.session(1).rec(4).blocks        = 7:8;
      subjectdata.session(1).rec(5).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:4;
      subjectdata.session(2).rec(3).blocks        = 5:6;
      subjectdata.session(2).rec(4).blocks        = 7:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 58

      subjectdata.name                    = 'BK';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:2;
      subjectdata.session(1).rec(3).blocks        = 3:5;
      subjectdata.session(1).rec(4).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = [2 3 4 6 7];
      subjectdata.session(2).rec(2).blocks        = 1:2;
      subjectdata.session(2).rec(2).block(2).missingtrials = 29:60;
      subjectdata.session(2).rec(3).blocks        = 2;
      subjectdata.session(2).rec(3).block(2).missingtrials = 1:45;
      subjectdata.session(2).rec(4).blocks        = 3:5;
      subjectdata.session(2).rec(6).blocks        = 6:8;
      subjectdata.session(2).rec(7).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;
      subjectdata.session(2).missingtrials        = 'rec2, stopped at trial 47. finished afterwards';

    case 59

      subjectdata.name                    = 'RO';
      subjectdata.drug                    = 'placebo';

      subjectdata.session(1).recsorder            = 2:7;
      subjectdata.session(1).rec(2).blocks        = 1;
      subjectdata.session(1).rec(3).blocks        = 2:3;
      subjectdata.session(1).rec(4).blocks        = 4:5;
      subjectdata.session(1).rec(5).blocks        = 6:7;
      subjectdata.session(1).rec(6).blocks        = 8:9;
      subjectdata.session(1).rec(7).blocks        = 10;
      subjectdata.session(1).rec(7).block(10).missingtrials = 1:50;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:2;
      subjectdata.session(2).rec(3).blocks        = 3:5;
      subjectdata.session(2).rec(4).blocks        = 6:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 60

      subjectdata.name                    = 'DL';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 3:4;
      subjectdata.session(1).rec(3).blocks        = 1:5;
      subjectdata.session(1).rec(4).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:4;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:7;
      subjectdata.session(2).rec(4).blocks        = 8:10;
      subjectdata.session(2).restingstate         = 1;

    case 61

      subjectdata.name                    = 'RW';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:3;
      subjectdata.session(1).rec(2).blocks        = 1:5;
      subjectdata.session(1).rec(3).blocks        = 6:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 62

      subjectdata.name                    = 'FZ';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:7;
      subjectdata.session(1).rec(2).blocks        = 1;
      subjectdata.session(1).rec(3).blocks        = 2:4;
      subjectdata.session(1).rec(4).blocks        = 5:6;
      subjectdata.session(1).rec(5).blocks        = 7:8;
      subjectdata.session(1).rec(6).blocks        = 9;
      subjectdata.session(1).rec(7).blocks        = 10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:5;
      subjectdata.session(2).rec(2).blocks        = 1:3;
      subjectdata.session(2).rec(3).blocks        = 4:6;
      subjectdata.session(2).rec(4).blocks        = 7:8;
      subjectdata.session(2).rec(5).blocks        = 9:10;
      subjectdata.session(2).restingstate         = 1;

    case 63

      subjectdata.name                    = 'HJ';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(2).block(3).missingtrials = 28:60;
      subjectdata.session(1).rec(3).blocks        = 3:6;
      subjectdata.session(1).rec(3).block(3).missingtrials = 1:27;
      subjectdata.session(1).rec(4).blocks        = 7:10;
      subjectdata.session(1).restingstate         = 1;
      subjectdata.session(1).missingtrials        = 'stopped block 3 at trial 27, continued from there in new recording';

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    case 64

      subjectdata.name                    = 'SE';
      subjectdata.drug                    = 'donepezil';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:4;
      subjectdata.session(1).rec(3).blocks        = 5:7;
      subjectdata.session(1).rec(4).blocks        = 8:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:6;
      subjectdata.session(2).rec(3).blocks        = 7:10;
      subjectdata.session(2).restingstate         = 1;

    case 65

      subjectdata.name                    = 'UR';
      subjectdata.drug                    = 'atomoxetine';

      subjectdata.session(1).recsorder            = 2:4;
      subjectdata.session(1).rec(2).blocks        = 1:3;
      subjectdata.session(1).rec(3).blocks        = 4:8;
      subjectdata.session(1).rec(4).blocks        = 9:10;
      subjectdata.session(1).restingstate         = 1;

      subjectdata.session(2).recsorder            = 2:3;
      subjectdata.session(2).rec(2).blocks        = 1:5;
      subjectdata.session(2).rec(3).blocks        = 6:10;
      subjectdata.session(2).restingstate         = 1;

    end

    % ====================================================================================== %
    % for each recording, get the name of the corresponding dataset
    % ====================================================================================== %

    if isnumeric(sj),
      allfiles = dir(subjectdata.rawdir);
      for session = 1:length(subjectdata.session),
        for rec = subjectdata.session(session).recsorder,
          expression = ['P\d*-S\d*' num2str(session) '_DotsPL_\d*_0' num2str(rec)];
          try
            name = allfiles(cellfun(@(x) ~isempty(regexp(x, expression, 'ignorecase')), {allfiles.name})).name;
            subjectdata.session(session).rec(rec).dataset = name;
          end
        end
      end
    end

    % =========================================== %
    % DATA FROM TABLE, copied from Excel
    % done by Eline
    % =========================================== %
    try
      t = readtable('tableSubjectData.csv');
      assert(strcmp(subjectdata.drug, t.drug(find(t.sj == sj))), 'drug does not match');
      flds = t.Properties.VariableNames;
      for f = 3:length(flds)-4,
        subjectdata.(flds{f})      = t.(flds{f})(find(t.sj == sj));
      end

      % heartrate, entered by Sven
      for session = 1:length(subjectdata.session),
        subjectdata.session(session).heartrate = [t.(['Heartrate_S' num2str(session) '_1'])(find(t.sj == sj))...
        t.(['Heartrate_S' num2str(session) '_2'])(find(t.sj == sj))];
      end
    end

  end
end
