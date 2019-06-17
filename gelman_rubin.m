function gelman_rubin

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

global mypath datasets datasetnames 

% GOODNESS OF FIT OF HIERARCHICAL MODEL
allgr = [];
for d = 1:length(datasets),

    gr = readtable(sprintf('%s/%s/stimcoding_nohist/gelman_rubin.txt', mypath, datasets{d}));
    gr = gr(~contains(gr.Var1, 'subj'), :)

    % grab group-level estimates
    allgr = [allgr; gr{:, 2}];
    if any(gr{:, 2} > 2),
    	disp(datasetnames{d})
    end
end

fprintf('Gelman-Rubin R-hat statistics: min %f, max %f, mean %f, std %f', ...
	min(allgr), max(allgr), mean(allgr), std(allgr));

% GOODNESS OF FIT OF G-SQ NONHIERARCHICAL MODEL

end