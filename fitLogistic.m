function [bias, slope, lapseLow, lapseHigh] = fitLogistic(x,y)

% make gamma and lambda symmetrical
pBest = fminsearchbnd(@(p) logistic_LL(p, ...
    x, y), [0 1 0.1 0.1], [-6 0 0 0], [6 20 1 1]);

bias        = pBest(1);
slope       = pBest(2);
lapseLow    = pBest(3);
lapseHigh   = pBest(3);

end

function err = logistic_LL(p, intensity, responses)
% see http://courses.washington.edu/matlab1/Lesson_5.html#1

% compute the vector of responses for each level of intensity
w   = logistic(p, intensity);

% negative loglikelihood, to be minimised
err = -sum(responses .*log(w) + (1-responses).*log(1-w));

end