% Simulates behaviour from timescale adaptation task using two independent
% leaky accumulators, one for each each alternative
%
% Inputs:
%       params(1) = beta (brightness exponent)
%       params(2) = sigmaA (accumulator noise)
%       params(3) = theta (decision bound)
%       params(4) = leak
%       params(5) = Ter (non-decision time)
%       
%       s = signal strength
%       d = signal duration
%       nsims = number of simulations to run

function ACC = sim2acc(params,s,d,nsims)
% Fixed parameters
dt = 1/60;             % model time-step (s)

baseL = 0.35;          % baseline stimulus input
sigmaL = 0.13;         % noise in stimulus input

trial_length = 5;
onset_range = [0.6 3.5];   % range of signal onsets (s)

% Free parameters
beta = params(1);       % brightness exponent (scales stimulus input)
sigmaA = params(2);     % 'internal' accumulator noise
theta = params(3);      % response threshold
leak = params(4);       % accumulator leakage
Ter = params(5);        % non-decision time

% Construct stimuli
input = (randn(2,nsims,(trial_length-Ter)/dt).*sigmaL)+baseL; % Raw noise-only input for each luminance patch
if s>0   % Adding signal to one luminance patch if required
    onsets = (onset_range(1)+(onset_range(2)-onset_range(1)).*rand(nsims,1))-Ter;  % drawing uniformly distributed signal onset latencies
    onsets(round(onsets./dt)==0) = dt;  % just in case onset was 0ms, which makes script error
    for i=1:nsims;
        input(1,i,round(onsets(i)/dt):(round(onsets(i)/dt)+(d/dt))) = input(1,i,round(onsets(i)/dt):(round(onsets(i)/dt)+(d/dt)))+s;
    end
end
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

ACC = zeros(1,nsims);  % prespecifying vector of accuracy markers
if s==0
    ACC(max(max(y,[],3))<theta) = 31;  % correct rejections
    ACC(max(max(y,[],3))>=theta) = 51;  % false alarms
else
    [I1,J1] = find(squeeze(y(1,:,:))>=theta);  % getting latencies of all passage times for each accumulator/trial
    [I2,J2] = find(squeeze(y(2,:,:))>=theta);
    
    [I1,Iindex] = unique(I1,'first'); J1 = J1(Iindex);  % keeping only *first* passage times for each accumulator/trial
    [I2,Iindex] = unique(I2,'first'); J2 = J2(Iindex);
    
    J1 = (J1.*dt)-onsets(I1);  % convert first passage times to RTs relative to signal onsets
    RT1 = nan(1,nsims); RT1(I1) = J1;  % and plug into vector with length = nsims
    
    J2 = (J2.*dt)-onsets(I2);
    RT2 = nan(1,nsims); RT2(I2) = J2;
    
    ts = 1:nsims;  % vector of trial numbers for later indexing
    dual_pass = find(ismember(ts,I1) & ismember(ts,I2));  % vector of trial numbers with dual bound crossings for later indexing
    %ACC(~ismember(ts,I1) & ~ismember(ts,I2)) = 41;  % miss
    %ACC((ismember(ts,I1) & RT1<0) | (ismember(ts,I2) & RT2<0)) = 71;  % premature responses
    %ACC((ismember(ts,I1) & RT1>d & ~ismember(ts,I2)) | (ismember(ts,I2) & RT2>d & ~ismember(ts,I1)) | (ismember(ts,dual_pass) & RT1>d & RT2>d)) = 81;  % too slow responses
    %ACC((ismember(ts,I2) & RT2>=0 & RT2<=d) & (~ismember(ts,I1) | ismember(ts,dual_pass(RT2(dual_pass)-RT1(dual_pass)<0)))) = 91;  % mislocalization
    ACC((ismember(ts,I1) & RT1>=0 & RT1<=d) & (~ismember(ts,I2) | ismember(ts,dual_pass(RT1(dual_pass)-RT2(dual_pass)<=0)))) = 1;  % hit                
end
