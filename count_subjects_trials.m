function count_subjects_trials

global mypath datasets

for d = 1:length(datasets)-1
    
    filename = dir(sprintf('%s/%s/*.csv', mypath, datasets{d}));
    dat  = readtable(sprintf('%s/%s/%s', mypath, datasets{d}, filename.name));
    
    nsubj(d) = numel(unique(dat.subj_idx));
    ntr(d) = height(dat);
end

fprintf('Total %d subjects, %d trials', sum(nsubj), sum(ntr));


% for MEG
d = 4;
dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
dat = dat(dat.session == 0, :);
for sj = unique(dat.subjnr)', sjdat = subjectspecifics(sj); gender(sj) = sjdat.gender; end

womencnt = 0; mencnt = 0;
for i = 1:length(gender),
    disp(gender{i});
    if ~isempty(gender{i}),
        if contains(gender{i}, 'F'),
            womencnt = womencnt + 1;
        elseif contains(gender{i}, 'M'),
            mencnt = mencnt + 1;
        end
    end
end
fprintf('Total %d men, %d women', mencnt, womencnt);
