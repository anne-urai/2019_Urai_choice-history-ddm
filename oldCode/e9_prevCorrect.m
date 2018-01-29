function e9_prevCorrect
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));
global mypath datasets datasetnames

close all;
for d = 1:length(datasets),
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    disp(datasets{d});
    
    switch d
        case {1 3 7} % only 1 coherence level
            results.z_prevcorrect = ...
                invlogit( results.z_1_1__stimcodingdczprevcorrect) - invlogit(results.z_1_2__stimcodingdczprevcorrect);
            results.z_preverror = ...
                invlogit( results.z_0_1__stimcodingdczprevcorrect) - invlogit(results.z_0_2__stimcodingdczprevcorrect);
            
            results.v_prevcorrect = ...
                results.dc_1_1__stimcodingdczprevcorrect - results.dc_1_2__stimcodingdczprevcorrect;
            results.v_preverror= ...
                results.dc_0_1__stimcodingdczprevcorrect - results.dc_0_2__stimcodingdczprevcorrect;
            
        case {2 4 5 6} % something weird happened with the coding of the 1-1, redo
            results.z_prevcorrect = ...
                invlogit( results.z_1_1__stimcodingdczprevcorrect) - invlogit(results.z_1_2__stimcodingdczprevcorrect);
            results.z_preverror = ...
                invlogit( results.z_c10__stimcodingdczprevcorrect) - invlogit(results.z_0_2__stimcodingdczprevcorrect);
            
            results.v_prevcorrect = ...
                results.dc_1_1__stimcodingdczprevcorrect - results.dc_1_2__stimcodingdczprevcorrect;
            % something weird
            results.v_preverror= ...
                results.dc_c10__stimcodingdczprevcorrect - results.dc_0_2__stimcodingdczprevcorrect;
    end
    
    % make nice-looking scatter plots
    bestcolor = linspecer(4, 'qualitative');
    colors = linspecer(5); % red blue green
    
    close all;
    subplot(4,4,1); hold on;
    scatterHistDiff(results.v_preverror, results.v_prevcorrect, [], [], colors(4, :));
    title(datasetnames{d});
    xlabel('v_{bias} after error');
    
    if ismember(d, [1 5]),
        ylabel('v_{bias} after correct');
    end
    
    subplot(4,4,5); hold on;
    scatterHistDiff(results.z_preverror, results.z_prevcorrect, [], [], colors(5, :));
    xlabel('z_{bias} after error');
    
    if ismember(d, [1 5]),
        ylabel('z_{bias} after correct');
    end
    tightfig;
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/prevCorrect_d%d.pdf', d));
    
end
close all;

end
