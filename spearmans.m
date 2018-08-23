

function [r, ci, pval, bf] = spearmans(a,b)

%[r1,pval1] = corr(a,b, 'type', 'spearman');

[r,pval,ci1,ci2] = corrcoef2(a,b,'spearman');
ci = [ci1 ci2];

% bayes factor of the correlation coefficient
bf = corrbf(r, numel(a));

end