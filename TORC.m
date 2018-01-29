%% read in a TORC audiofile from Jan-Willem and plot the signal and noise timecourses

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

cd ~/Downloads/sounds/

[target, Fs]      = audioread('TORC_TARGET.wav');
[noise, Fs2]      = audioread('TORC_424_01_h501.wav');
assert(Fs == Fs2, 'sampling rate is different between signal and noise');
targetnoise = target(1:length(noise)) + noise;

close all;
subplot(211);
plot(targetnoise, 'k');
axis tight; xlim([1 600]);
axis off;

subplot(212);
plot(noise, 'k');
axis tight; xlim([1 600]);
axis off;
print(gcf, '-dpdf', 'TORC.pdf');

% also print the spectograms, as in McGinley et al Neuron
subplot(441);
spectrogram(targetnoise, 400, 300, [], Fs, 'yaxis');
xlabel('Time (s)');
ylabel('Freq (kHz)');
