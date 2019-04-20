function bias_by_congruency

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;

for d = 1:length(datasets),


    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);


    % {'z_0_2__stimcodingdczprevrespcongruency'     }
    % {'z_0_1__stimcodingdczprevrespcongruency'     }
    % {'z_1_2__stimcodingdczprevrespcongruency'     }
    % {'z_1_1__stimcodingdczprevrespcongruency'     }
    % {'dc_0_2__stimcodingdczprevrespcongruency'    }
    % {'dc_0_1__stimcodingdczprevrespcongruency'    }
    % {'dc_1_2__stimcodingdczprevrespcongruency'    }
    % {'dc_1_1__stimcodingdczprevrespcongruency'    }

    %% STARTING POINT
    z_congruent = results.z_1_2__stimcodingdczprevrespcongruency - results.z_1_1__stimcodingdczprevrespcongruency;
    z_incongruent = results.z_0_2__stimcodingdczprevrespcongruency - results.z_0_1__stimcodingdczprevrespcongruency;

    close all;
    sp = subplot(4,4,1);
    hold on;

    plot([z_congruent z_incongruent]');
    pval = permtest(z_congruent, z_incongruent);
    mysigstar(gca, [1 2], max(get(gca, 'ylim')), pval);
    legtext{d} = cat(2, datasetnames{d}{1}, ' ', datasetnames{d}{2});
    set(gca, 'xtick', [1 2], 'xticklabel', {'congruent', 'incongruent'})
    title(datasetnames{d});

    ylabel({'Starting point shift'});
    offsetAxes;
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/congruency_z_%d.pdf', d));


    %% STARTING POINT
    dc_congruent = results.dc_1_2__stimcodingdczprevrespcongruency - results.dc_1_1__stimcodingdczprevrespcongruency;
    dc_incongruent = results.dc_0_2__stimcodingdczprevrespcongruency - results.dc_0_1__stimcodingdczprevrespcongruency;

    close all;
    sp = subplot(4,4,1);
    hold on;

    plot([dc_congruent dc_incongruent]');
    pval = permtest(dc_congruent, dc_incongruent);
    mysigstar(gca, [1 2], max(get(gca, 'ylim')), pval);
    title(datasetnames{d});
    set(gca, 'xtick', [1 2], 'xticklabel', {'congruent', 'incongruent'})

    ylabel({'Drift bias shift'});
    offsetAxes;
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/congruency_vbias_%d.pdf', d));
end
