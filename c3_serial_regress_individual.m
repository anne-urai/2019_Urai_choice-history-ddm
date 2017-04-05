function c3_serial_regress_individual
% using the stimcoding models where dc and/or z depend on previous
% response, plot

clear;
datasets{1} = '/Users/anne/Data/projects/0/neurodec/Data/MEG-PL/HDDM/summary';
datasets{2} = '~/Data/RT_RDK/HDDM/summary';

for d = 1:2,
    
    usepath = datasets{d};
    close all;
    fz = 8;
    set(groot, 'defaultaxesfontsize', fz, ...
        'defaultaxestitlefontsizemultiplier', 1, 'defaultaxestitlefontweight', 'normal');
    
    % get the individual and group stuff
    dat_dc      = load(sprintf('%s/regress_dc_prevresp_prevpupil_prevrt_all.mat', usepath));
    dat_z       = load(sprintf('%s/regress_z_prevresp_prevpupil_prevrt_all.mat', usepath));
    % dat_both    = readtable(sprintf('%s/regress_dc_z_prevresp_group_traces_concat.csv', usepath));
    
    
    subplot(221);
    scatter(dat_dc.individuals.v_prevresp_mean, dat_z.individuals.z_prevresp_mean);
    title('Repetition'); lsline;
    
    subplot(222);
    scatter(dat_dc.individuals.v_prevpupil_prevresp_mean, dat_z.individuals.z_prevpupil_prevresp_mean);
    title('Pupil modulation'); lsline;
    
    subplot(223);
    scatter(dat_dc.individuals.v_prevrt_prevresp_mean, dat_z.individuals.z_prevrt_prevresp_mean);
    title('RT modulation'); lsline;
    
    suptitle(usepath);
    print(gcf, '-dpdf', sprintf('%s/serial_regress.pdf', usepath));
end

end
