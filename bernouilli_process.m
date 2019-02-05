function rep_trivial = bernouilli_process(lag, prep)

% define each trial as a repeat 1 or alternate 0 with probability p
% (for each individual, p will be computed as in Figure 2b)
%
% Q: given p, what is the probability of repeating the choice that was made 
% l lags into the past?
%

% example:
% 
% lags = 1:10;
% prep = 0.1:0.1:0.9;
% for l = 1:length(lags),
%     for r = 1:length(prep),
%         trivial(l,r) = bernouilli_process(lags(l), prep(r));
%     end
% end
% figure; subplot(221); imagesc(lags, prep,  trivial');
% xlabel('Lags'); ylabel('Generative P(repeat) at lag 1');
% c = colorbar; c.Label.String('Expected P(repeat)');

% ======================================================================== %
% first, enumerate all the possible sequences and their likelihood given p
% ======================================================================== %
%{
lag 1
XY = A = (1-p)
XX = R = p

lag 2
XXX = RR = p*p
XXY = RA = p*(1-p)
XYX = AA = (p-1)*(p-1)
XYY = AR = (p-1)*p

lag 3
XXXX = RRR = p*p*p
XXXY = RRA
XXYX = RAA
XXYY = RAR
XYXX = ARR
XYXY = AAA
XYYX = ARA
XYYY = ARR

%}

if ~exist('lag', 'var'), lag = 7; end
if ~exist('prep', 'var'), prep = 0.6; end

% DO THIS automagically
% generate all the possible sequences of repetitions and alternations
all_binary_sequences = @(x) dec2bin(2^x-1:-1:0)-'0';
sequences_rep_alt = all_binary_sequences(lag);

% COMPUTE THE PROBABILITY OF EACH SEQUENCE
sequences_prob = sequences_rep_alt;
sequences_prob(sequences_prob == 1) = prep;
sequences_prob(sequences_prob == 0) = 1 - prep;

prob = prod(sequences_prob, 2);
assert(roundn(sum(prob), -3) == 1, 'sequence probabilities must sum up to 1');

% FOR EACH LAG, CHECK WHICH SEQUENCES HAVE REPEATS OF THE CURRENT TRIAL
% code alternates as -1
sequences_rep_alt(sequences_rep_alt == 0) = -1;
% take the product: two alternations cancel each other out
repeat = prod(sequences_rep_alt, 2);
repeat(repeat == -1) = 0;

% compute the expected repetition probability at this lag
rep_trivial = sum(prob(repeat == 1));


