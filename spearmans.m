

function [r, ci, pval, bf] = spearmans(a,b)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

[r,pval,ci1,ci2] = corrcoef2(a,b,'spearman');
ci = [ci1 ci2];

% bayes factor of the correlation coefficient
bf = corrbf(r, numel(a));

end