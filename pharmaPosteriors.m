function pharmaPosteriors

global mypath

% neutral vs biased plots
datasets = {'MEG_MEGsessions'};
close all;

for d = 1:length(datasets),
    
    traces = readtable(sprintf('%s/%s/stimcoding_dc_z_prevresp_pharma/all_traces.csv', mypath, datasets{d}));
    
    z_placebo = invlogit(traces.z_trans_placebo_1_) - invlogit(traces.z_trans_placebo__1_);
    z_atomox  = invlogit(traces.z_trans_atomoxetine_1_) - invlogit(traces.z_trans_atomoxetine__1_);
    z_donepezil  = invlogit(traces.z_trans_donepezil_1_) - invlogit(traces.z_trans_donepezil__1_);
    
    v_placebo = (traces.dc_placebo_1_) - (traces.dc_placebo__1_);
    v_atomox  = (traces.dc_atomoxetine_1_) - (traces.dc_atomoxetine__1_);
    v_donepezil  = (traces.dc_donepezil_1_) - (traces.dc_donepezil__1_);
    
    colors = [0 0 0; 1 0 0; 0 0 1];
    cb = cbrewer('qual', 'Set1', 4);
    colors = [0 0 0; cb(1:2, :)];
    subplot(441); hold on;
    histogram_smooth(z_placebo, colors(1, :));
    histogram_smooth(z_atomox, colors(2, :));
    histogram_smooth(z_donepezil, colors(3, :));
    xlabel('History shift in z');
    ylabel('Posterior probability');
    set(gca, 'yticklabel', []);
    box off;
    
    % ADD P-VALUES
    posteriorpval = @(dat1, dat2) min([mean(dat2 > dat1) mean(dat2 < dat1)]);
    pvalat = posteriorpval(z_atomox, z_placebo)
    pvalat = posteriorpval(z_donepezil, z_placebo)
    pvalat = posteriorpval(z_donepezil, z_atomox)
    
    subplot(442); hold on;
    h1 = histogram_smooth(v_placebo, colors(1, :));
    h2 = histogram_smooth(v_atomox, colors(2, :));
    h3 = histogram_smooth(v_donepezil, colors(3, :));
    xlabel('History shift in v');
    ylabel('Posterior probability');
    set(gca, 'yticklabel', []);
    box off;
    
    pvalat = posteriorpval(v_atomox, v_placebo)
    pvalat = posteriorpval(v_donepezil, v_placebo)
    pvalat = posteriorpval(v_donepezil, v_atomox)
    
    l = legend([h1 h2 h3], {'Placebo', 'Atomoxetine', 'Donepezil'});
    l.Box = 'off';
    l.Position(1) = l.Position(1) + 0.15;
    
    print(gcf, '-deps', sprintf('~/Data/serialHDDM/pharmaPosteriors.eps'));
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/pharmaPosteriors.pdf'));
    
end

end


function h = histogram_smooth(x, color2)

[f,xi] = ksdensity(x);
a1 = area(xi, f, 'edgecolor', 'none', 'facecolor', color2, 'facealpha', 0.4);

% % Make area transparent
% drawnow; % pause(0.05);  % This needs to be done for transparency to work
% a1.Face.ColorType = 'truecoloralpha';
% a1.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3

% area
h = plot(xi, f, 'color', color2, 'linewidth', 1);
set(gca, 'color', 'none');

end
