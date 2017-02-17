clear; close all; clc;
addpath(genpath('~/Code/Elisa'));

% get all those subjects who have a folder
cd ~/ElisaBattistoni;
s = dir('P*');
s = {s(:).name};
for i = 1:length(s), subjects(i) = str2num(s{i}(2:3)); end

for sj = 4; %unique(subjects),
    
    tic
    clearvars -except subjects sj

    subjectdata.files   = dir(sprintf('~/ElisaBattistoni/P%02d/behavioral/P%d_s*.mat', sj, sj));
    datetimeidx         = regexp({subjectdata.files(:).name}, '\d*', 'match');
    % these indices now contain sj nr, session nr, 
    % year, month, hour, minutes, seconds
    
    % do some quick checks on this
    for i = 1:length(datetimeidx),
        assert(length(datetimeidx{i}) == 8, 'filename does not have the right format');
        assert(strcmp(datetimeidx{i}(3), '2014')==1, 'could not find the right datestr in the filename');
        sortingidx(i) = str2double(cell2mat(datetimeidx{i}(3:end)));
    end
    
    % now take the indices and sort the filenames accordingly
    [~,idx] = sort(sortingidx);
    files = {subjectdata.files(idx).name}; % this returns the filenames sorted by the date created
    
    % alldata = struct('RT', [], 'correct', [], 'session', [], 'response', [], 'increment', []);
    alldata = [];
    
    % concatenate all data from this subject
    for f = 1:length(files),
        load([sprintf('~/ElisaBattistoni/P%02d/behavioral/', sj) files{f}]);
        % to save memory, clear the dots for now
        clearvars stim flip
        
        for b = 1:setup.nblocks,
            for t = 1:setup.ntrials,
                if ~isnan(results.response(b,t)),
                    % skip if there are no responses at all in this block

                    thistrial = [results.correct(b,t); results.response(b,t);
                        results.RT(b,t); dots.direction(b,t); setup.session; setup.threshold];
                    
                    switch isequal(results.response(b,t),dots.direction(b,t))
                        case 1
                            assert(results.correct(b,t) == 1, 'correctness coding error');
                        case 0
                            assert(results.correct(b,t) == 0, 'correctness coding error');
                    end
                            
                    
                    alldata = [alldata thistrial];
                    % alldata now contains: correct response RT
                    % direction session threshold
                end
            end
        end
    end
    
    %end
    
    disp(sprintf('Saving P%02d', sj));
    save(sprintf('~/ElisaBattistoni/P%02d/behavioral/P%d_alldata.mat', sj, sj), 'alldata');
end