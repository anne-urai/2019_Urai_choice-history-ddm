function [res, dat] = A5b_FruendModel(res)

res.choiceW = nan(size(res.subjnr));
res.stimW = nan(size(res.subjnr));
res.pupil_choiceW = nan(size(res.subjnr));
res.pupil_stimW = nan(size(res.subjnr));
res.rt_choiceW = nan(size(res.subjnr));
res.rt_stimW = nan(size(res.subjnr));

% ================================================================== %
% WRITE FILES
% ================================================================== %

addpath /Users/anne/Dropbox/code/MEG;
subjectdata = subjectspecifics('ga');
% this file should also include sessions without pupil measurements
alldata = readtable(sprintf('%s/Data/CSV/2ifc_pupildata_allsj.csv', subjectdata.path));

for s = 0:2, % all data, and session 1 and 2 separately
    for sj = unique(alldata.subjnr)',
        
        if exist(sprintf('%s/serial/txt/meg_pupil+rt_sj%02d_session%d.txt', subjectdata.path, sj, s), 'file');
            %   continue;
        end
        
        % use files with cleaner pupil data
        switch s
            case 0
                % both sessions together
                data   = alldata(alldata.subjnr == sj  & ~isnan(alldata.rt), :);
            otherwise
                data   = alldata(alldata.subjnr == sj & ~isnan(alldata.rt) & alldata.sessionnr == s , :);
        end
        if isempty(data), continue; end
        
        % generate block nrs, NOT identical to session nrs!
        % History effects should not continue beyond a block
        blockchange = find(diff(data.trialnr) < 0);
        blocknrs = zeros(height(data), 1);
        for b = 1:length(blockchange)-1,
            blocknrs(blockchange(b)+1:blockchange(b+1)) = blocknrs(blockchange(b))+1;
        end
        blocknrs(blockchange(end)+1:end) = blocknrs(blockchange(end))+1;
        
        % double modulation: decision pupil and rt
        if any(isnan(data.decision_pupil)),
            data.decision_pupil(isnan(data.decision_pupil)) = 0; % do not modulate
        end
        
        newdat = [blocknrs data.sessionnr abs(data.stim) (data.stim > 0) (data.response > 0) ...
            zscore(data.decision_pupil) zscore(log(data.rt + 0.2))];
        assert(~any(isnan(newdat(:))));
        dlmwrite(sprintf('%s/serial/txt/meg_pupil+rt_sj%02d_session%d.txt', subjectdata.path, sj, s), ...
            newdat,'delimiter','\t','precision',4);
    end
    
    % ================================================================== %
    % RUN THE FRUEND TOOLBOX
    % ================================================================== %
    
    mods = {'pupil+rt'};
    addpath(genpath('/Users/anne/Data/pupilUncertainty_FigShare/Code/'));
    addpath(genpath('/Users/anne/Data/pupilUncertainty_FigShare/Code/serial-dependencies/'));
    
    for m = 1:length(mods),
        
        if exist(sprintf('%s/serial/historyweights_%s_session%d.mat', subjectdata.path, mods{m}, s), 'file');
            continue;
        end
        
        cd('/Users/anne/Data/pupilUncertainty_FigShare/Code/serial-dependencies/');
        system([sprintf('for sj in {2..65}; do filename=$(printf "%s/serial/txt/meg_%s_sj%%02d_session%d.txt" $sj);', subjectdata.path, mods{m}, s), ...
            sprintf('echo $filename; python2.7 analysis.py -fr -n10 -p "%s/serial/outp/" $filename; sleep 5; done', subjectdata.path)]);
    end
    
    % get the output back into something Matlab can work with
    cd /Users/anne/Dropbox/code/MEG/Behaviour;
    retrieveDataFromPython(mods{m}, subjectdata.path, s);
    
    % ================================================================== %
    % ANALYZE THE DATA
    % ================================================================== %
    
    load(sprintf('%s/serial/historyweights_%s_session%d.mat', subjectdata.path, mods{m}, s));
    res.choiceW(res.session == s)       = dat.response(subjectdata.all, 1);
    res.stimW(res.session == s)         = dat.stimulus(subjectdata.all, 1);
    res.pupil_choiceW(res.session == s) = dat.response_pupil(subjectdata.all, 1);
    res.rt_choiceW(res.session == s)    = dat.response_rt(subjectdata.all, 1);
    res.pupil_stimW(res.session == s)   = dat.stimulus_pupil(subjectdata.all, 1);
    res.rt_stimW(res.session == s)      = dat.stimulus_rt(subjectdata.all, 1);
    
end
