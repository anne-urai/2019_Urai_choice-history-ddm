

%% GET ANKE'S DATA, AND DERIVE A COHERENCE BY DRIFT RATE FUNCTION
global datasets
d = 4;
load(sprintf('%s/summary/%s/stimcoding_nohist_all.mat', mypath, datasets{d}));

dat.coherence_levels = [0 3 9 27 81];
dat.drift_rates      = [individuals.v_c0_mean' ...
    individuals.v_c3_mean' ...
    individuals.v_c9_mean' ...
    individuals.v_c27_mean' ...
    individuals.v_c81_mean' ];

close;
plot(dat.coherence_levels, dat.drift_rates)
xlabel('% coherence');
ylabel('Drift rate');
print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/Anke_driftrates.pdf',d));


save(sprintf('%s/summary/%s/driftRates.mat', mypath, datasets{d}));
