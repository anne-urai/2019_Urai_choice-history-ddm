% Calculates likelihood term for given trial-type (noise-only,
% signal-present) given by 2-accumulator LCA model
%
% Inputs:
%      data(1:2) = [Nnoise Nfa] or [Nsignal Nhit], depending on trial-type
%      simrate = simulated FA/hit rate, depending on trial-type

function L = L2acc(data,simrate)

% [~,L] = evalc('nchoosek(data(1),data(2))*(simrate.^data(2))*((1-simrate).^(data(1)-data(2)));');

% Running initial calculations in log space and then taking the exponent -
% otherwise nchoosek runs into numerical precision issues
L = exp((gammaln(data(1)+1)-gammaln(data(1)-data(2)+1)-gammaln(data(2)+1))+log(simrate.^data(2))+log((1-simrate).^(data(1)-data(2))));
