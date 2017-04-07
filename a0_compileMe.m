function [] = compileMe(fname)

% make my life easier
if ~strcmp(fname(end-1:end), '.m'); fname = [fname '.m']; end % add extension
disp(fname);

% these paths will be added at compilation
addpath(genpath('~/code/Tools'));
addpath(genpath('~/code/Tools'));
addpath('~/Documents/fieldtrip/');
ft_defaults;

% options: compile verbose, only use the toolboxes we really need
% !!! runtime options should be preceded by - to work!
% dont need to activate the -nojvm flag, can still plot from executable
% no need to specify additional files
mcc('-mv', '-N', '-p', 'stats', '-p', 'signal', '-p', 'curvefit',...
    '-R', '-nodisplay', '-R', '-singleCompThread', ...
    fname);

delete mccExcludedFiles.log
delete run*.sh
delete readme.txt
delete requiredMCRProducts.txt

%% generate allparams file to add to st
subjects =  [3     4     5     6     7     8     9    10    11    12    13   ...
    14    15    17    18    19    20    21    22    23    24    25];
alldat = [];
for sj = subjects,
  for session = 1:5,
    for block = 1:10,
      alldat = [alldat; sj session block];
    end
  end
end

% write to a file
dlmwrite('allblocks', alldat, 'delimiter', ' ');

% ============================================ #
% parameter file for HDDM
% ============================================ #

datasets   = 4:6; % RT-RDK and MEG-PL
models     = [2 7 8]; % the nr of the models
nrTraces   = 30; % nr of chains, 15 cores/node (so make sure this is a multiple of 15)

alldat = [];
  for c = 0:nrTraces-1,
    for a = datasets,
      for b = models,
      alldat = [alldat; a b c];
    end
  end
end

% write to a file
dlmwrite('hddmparams', alldat, 'delimiter', ' ');
size(alldat)
