
function [] = a0_parameterFiles()

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
datasets   = [0:5]; % RT-RDK and MEG-PL
models     = [9]; % the nr of the models
nrTraces   = 30; % nr of chains, 15 cores/node (so make sure this is a multiple of 15)

alldat = [];
for n = nsmp,
    for b = models,
        for a = datasets,
            for c = 0:nrTraces-1, % put all chains of same model together on a node
				if ~(a == 4 && b == 4),
                alldat = [alldat; a b c n];
				end
            end
        end
    end
end

%% add the MEG st, will take forever so need more parallel chains
alldat = [];
for c = 0:29,
	alldat = [alldat; 1 8 c nsmp];
end
for c = 0:60-1,
alldat = [alldat; 4 4 c 1000];
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
