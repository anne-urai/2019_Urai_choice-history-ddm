
function [] = make_stopos_file()

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

% ============================================ #
% parameter file for HDDM
% ============================================ #s

nsmp       = [5000]
datasets   = [4]; % dataset number
models     = [11]; % the nr of the models
nrTraces   = 30; % nr of chains, 15 cores/node (so make sure this is a multiple of 15)

alldat = [];
for n = nsmp,
    for b = models,
        for a = datasets,
            for c = 0:nrTraces-1, % put all chains of same model together on a node
                alldat = [alldat; a b c n];
            end
        end
    end
end

% write to a file
dlmwrite('hddmparams', alldat, 'delimiter', ' ');
size(alldat)

% PPC
alldat = [];
for a = datasets,
    for v = models,
        alldat = [alldat; a v];
    end
end

% write to a file
dlmwrite('hddmparams_PPC', alldat, 'delimiter', ' ');
size(alldat)
%% for PPC & chiSquare
