% Simulates behaviour from timescale adaptation task using two independent
% leaky accumulators, one for each each alternative
%
% Inputs:
%       params(1) = beta (brightness exponent)
%       params(2) = sigmaA (accumulator noise)
%       params(3) = leak
%       params(4) = Ter (non-decision time)
%       
%       nsims = number of simulations to run

function y = sim2acc_no_bound(params,nsims)

% Fixed parameters
dt = 1/60;             % model time-step (s)

baseL = 0.35;          % baseline stimulus input
sigmaL = 0.11;         % noise in stimulus input

trial_length = 5;

% Free parameters
beta = params(1);       % brightness exponent (scales stimulus input)
sigmaA = params(2);     % 'internal' accumulator noise
leak = params(3);       % accumulator leakage
Ter = params(4);        % non-decision time

% Construct stimuli
input = (randn(2,nsims,(trial_length-Ter)/dt).*sigmaL)+baseL; % Raw noise-only input for each luminance patch
input(input<0) = 0;  % setting any negative luminance values to zero
input = (input.^beta)-(baseL.^beta); % Transforming input using brightness exponent and subtracting baseline luminance

% Simulating behaviour
noise = randn(2,nsims,(trial_length-Ter)/dt).*sigmaA;  % pre-drawing noise to speed things up
y = zeros(2,nsims,1+ceil((trial_length-Ter)/dt));
for t = 1:ceil((trial_length-Ter)/dt)  % looping through time only, simulating all trials simultaneously
    y(1,:,t+1) = (y(1,:,t).*(1-leak))+input(1,:,t)+noise(1,:,t);
    y(2,:,t+1) = (y(2,:,t).*(1-leak))+input(2,:,t)+noise(2,:,t);
end
y = y(:,:,2:end);  % dumping t=0