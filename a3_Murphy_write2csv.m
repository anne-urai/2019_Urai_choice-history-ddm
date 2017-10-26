mypath = '~/Data/PeterMurphy/MATs';
cd(mypath);
clear;

subjects    = dir('*_PMF.mat');
subjects    = regexp({subjects(:).name}, '\d*', 'match');
subjects    = cellfun(@str2num, [subjects{:}]);

%% read in all data
for sj = subjects
    for session = 1:5
        load(sprintf('%d_%d.mat', sj, session));
        
        thistab = array2table(data_block(:, [1 2 4 5 7]), 'variablenames', ...
            {'coherence', 'correct', 'rt', 'stimulus', 'isi'});
        thistab.subj_idx    = sj * ones(size(thistab.rt));
        thistab.session     = session * ones(size(thistab.rt));
        thistab.trial       = transpose(1:100);
        
        if sj == subjects(1) && session == 1
            t = thistab;
        else
            t = cat(1, t, thistab);
        end
        
    end
end

%% recode
t.stimulus = sign(t.stimulus - 10);
t.response = t.stimulus;
t.response(t.correct == 0) = - t.stimulus(t.correct == 0);
t.rt(t.rt < 0.2) = NaN;

% test
assert(isequal(t.correct, (t.stimulus == t.response)));
t.prevresp = circshift(t.response, 1);
t.prev2resp = circshift(t.response, 2);
t.prev3resp = circshift(t.response, 3);

t.prevstim = circshift(t.stimulus, 1);
t.prev2stim = circshift(t.stimulus, 2);
t.prev3stim = circshift(t.stimulus, 3);

t.prevrt   = circshift(nanzscore(log(t.rt)), 1);

% remove trials without previous
% remove trials where the previous trial was not immediately preceding
wrongtrls               = find([NaN; diff(t.trial)] ~= 1);
t(wrongtrls, :)         = [];
t(isnan(t.rt), :)       = [];
t.response(t.response == -1) = 0;
t.coherence             = [];

writetable(t, sprintf('~/Data/HDDM/Murphy/murphy_hddm.csv'));
disp('done!');

      