function kernels_lags_regression

global mypath datasets datasetnames
addpath(genpath('~/code/Tools'));
warning off; close all;
colors = cbrewer('qual', 'Set2', length(datasets));

%% PREVIOUS STIMULUS AND RESPONSE

mdls = {'regress_z_lag3', 'regress_dc_lag3', 'regress_dcz_lag3'};
types = {'stim', 'resp'};

for m = 1:length(mdls),
    for t = 1:length(types),
        
        close all;
        switch m
            case 1
                
                subplot(4,3,1); hold on;
                for d = 1:length(datasets),
                    load(sprintf('%s/summary/%s/%s_all.mat', ...
                        mypath, datasets{d}, mdls{m}));
                    mat_dc = [group.(sprintf('z_prev%s_mean', types{t})) ...
                        group.(sprintf('z_prev2%s_mean', types{t})) ...
                        group.(sprintf('z_prev3%s_mean', types{t}))];
                    plot(1:3, mat_dc, 'color', colors(d, :));
                end
                switch types{t}
                    case 'stim'
                        ylabel('z ~ previous stimulus');
                    case 'resp'
                        ylabel('z ~ previous respones');
                        
                end
                xlabel('Lag (past trials)');
                
                drawnow; tightfig;
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_m%d_type%d.pdf', m, t));
                
            case 2
                
                subplot(4,3,1); hold on;
                for d = 1:length(datasets),
                    load(sprintf('%s/summary/%s/%s_all.mat', ...
                        mypath, datasets{d}, mdls{m}));
                    mat_dc = [group.(sprintf('dc_prev%s_mean', types{t})) ...
                        group.(sprintf('dc_prev2%s_mean', types{t})) ...
                        group.(sprintf('dc_prev3%s_mean', types{t}))];
                    plot(1:3, mat_dc, 'color', colors(d, :));
                end
                switch types{t}
                    case 'stim'
                        ylabel('dc ~ previous stimulus');
                    case 'resp'
                        ylabel('dc ~ previous respones');
                        
                end
                xlabel('Lag (past trials)');
                
                drawnow; tightfig;
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_m%d_type%d.pdf', m, t));
                
                
                
            case 3
                
                subplot(4,3,1); hold on;
                for d = 1:length(datasets),
                    load(sprintf('%s/summary/%s/%s_all.mat', ...
                        mypath, datasets{d}, mdls{m}));
                    mat_dc = [group.(sprintf('dc_prev%s_mean', types{t})) ...
                        group.(sprintf('dc_prev2%s_mean', types{t})) ...
                        group.(sprintf('dc_prev3%s_mean', types{t}))];
                    plot(1:3, mat_dc, 'color', colors(d, :));
                end
                switch types{t}
                    case 'stim'
                        ylabel('dc ~ previous stimulus');
                    case 'resp'
                        ylabel('dc ~ previous respones');
                        
                end
                xlabel('Lag (past trials)');
                
                drawnow; tightfig;
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_m%d_dc_type%d.pdf', m, t));
                
                
                % shared model, x
                close all; subplot(4,3,1); hold on;
                for d = 1:length(datasets),
                    load(sprintf('%s/summary/%s/%s_all.mat', ...
                        mypath, datasets{d}, mdls{m}));
                    mat_dc = [group.(sprintf('z_prev%s_mean', types{t})) ...
                        group.(sprintf('z_prev2%s_mean', types{t})) ...
                        group.(sprintf('z_prev3%s_mean', types{t}))];
                    plot(1:3, mat_dc, 'color', colors(d, :));
                end
                switch types{t}
                    case 'stim'
                        ylabel('z ~ previous stimulus');
                    case 'resp'
                        ylabel('z ~ previous respones');
                        
                end
                xlabel('Lag (past trials)');
                
                drawnow; tightfig;
                print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/regressionkernels_m%d_z_type%d.pdf', m, t));
                
        end
        
    end
end