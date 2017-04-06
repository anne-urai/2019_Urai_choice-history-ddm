
function serial_regress
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
    
    % first, group-level traces
    dat_dc      = readtable(sprintf('%s/regress_dc_prevresp_group_traces_concat.csv', usepath));
    dat_z       = readtable(sprintf('%s/regress_z_prevresp_group_traces_concat.csv', usepath));
    %dat_both    = readtable(sprintf('%s/regress_dc_z_prevresp_group_traces_concat.csv', usepath));
    
    % top row, separate dc and z
    sp1 = subplot(441);
    plothist(dat_dc, 'v_prevresp'); %xlim([-0.05 0.08])
    
    sp2 = subplot(443); hold on;
    plothist(dat_z, 'z_prevresp'); %xlim([-0.04 0.04]);
    
    %% now also the prevresp and prevpupil traces
    % first, group-level traces
    dat_dc      = readtable(sprintf('%s/regress_dc_prevresp_prevpupil_prevrt_group_traces_concat.csv', usepath));
    dat_z       = readtable(sprintf('%s/regress_z_prevresp_prevpupil_prevrt_group_traces_concat.csv', usepath));
    
    % top row, separate dc and z
    subplot(445);
    plothist(dat_dc, 'v_prevresp'); %xlim([-0.05 0.08])
    set(gca, 'xlim', get(sp1, 'xlim'));
    
    subplot(447); hold on;
    plothist(dat_z, 'z_prevresp');  %xlim([-0.04 0.04]);
    set(gca, 'xlim', get(sp2, 'xlim'));
    
    % now the modulatory terms themselves
    subplot(4,4,9);
    plothist(dat_dc, 'v_prevpupil_prevresp');
    
    subplot(4,4,10); hold on;
    plothist(dat_dc, 'v_prevrt_prevresp');
    
    subplot(4,4,11); hold on;
    plothist(dat_z, 'z_prevpupil_prevresp');
    
    subplot(4,4,12); hold on;
    plothist(dat_z, 'z_prevrt_prevresp');
    
    suptitle(usepath);
    print(gcf, '-dpdf', sprintf('%s/serial_regress.pdf', usepath));
end

end

function plothist(dat, param)

histogram(dat.(param), 'displaystyle', 'stairs');
axis tight; box off;
xlabel(param, 'interpreter', 'none'); %xlim([-0.2 0.2]);
hold on;
plot([0 0], get(gca, 'ylim'), 'k');

end

