
% using the stimcoding models where dc and/or z depend on previous
% response, plot

usepath = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/HDDM/summary';
close all;
fz = 8;
set(groot, 'defaultaxesfontsize', fz, ...
    'defaultaxestitlefontsizemultiplier', 1, 'defaultaxestitlefontweight', 'normal');

% first, group-level traces
dat_dc      = readtable(sprintf('%s/stimcoding_prevresp_dc_group_traces_concat.csv', usepath));
dat_z       = readtable(sprintf('%s/stimcoding_prevresp_z_group_traces_concat.csv', usepath));
dat_both    = readtable(sprintf('%s/stimcoding_prevresp_dc_z_group_traces_concat.csv', usepath));

% top row, separate dc and z
subplot(221); hold on;
histogram(dat_dc.dc_1, 'displaystyle', 'stairs');
histogram(dat_dc.dc_2, 'displaystyle', 'stairs');
axis tight; box off; 
xlabel('dc'); xlim([-0.2 0.2]);
ylabel('Separate models');

subplot(222); hold on;
histogram(dat_z.z_1, 'displaystyle', 'stairs');
histogram(dat_z.z_2, 'displaystyle', 'stairs');
axis tight; box off; 
xlabel('z'); xlim([0.45 0.51]);

% top row, separate dc and z
subplot(223); hold on;
histogram(dat_both.dc_1, 'displaystyle', 'stairs');
histogram(dat_both.dc_2, 'displaystyle', 'stairs');
axis tight; box off; 
xlabel('dc'); xlim([-0.2 0.2]);
ylabel('Joint model');

subplot(224); hold on;
histogram(dat_both.z_1, 'displaystyle', 'stairs');
histogram(dat_both.z_2, 'displaystyle', 'stairs');
axis tight; box off; 
xlabel('z'); xlim([0.45 0.51]);

print(gcf, '-dpdf', sprintf('%s/serial_stimcoding.pdf', usepath));


