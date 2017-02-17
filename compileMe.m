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
    14    15    17    18    19    20    21    22    23    24    25]';

% write to a file
dlmwrite('allsubjects', subjects, 'delimiter', ' ');
