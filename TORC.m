cd ~/Downloads/sounds/

target = audioread('TORC_TARGET.wav');
noise = audioread('TORC_424_01_h501.wav');

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
