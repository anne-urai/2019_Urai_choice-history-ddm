function test_regression_recoding

global mypath datasets
addpath(genpath('~/code/Tools'));
warning off; close all;

for d = 1:length(datasets),
    dat = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    dat = dat(dat.session == 0, :);
    %dat = readtable('auditory_yesno_hddmfits.csv');
    
    % get the true correct/error weights for lag 7
    vars = {'z_correct', 'z_error', 'v_correct', 'v_error', ...
        'z_error_recode', 'z_correct_recode', 'v_correct_recode', 'v_error_recode'};
    numlags = 7;
    for m = 1:length(vars),
        alldata.(vars{m})   = nan(length(unique(dat.subjnr)), numlags);
        alldata.([vars{m} '_pval'])   = nan(length(unique(dat.subjnr)), numlags);
    end
    bestmodelname = 'regressdczlag7';
    
    for l = 1:numlags,
        if l == 1,
            lname = [];
        else
            lname = l;
        end
        
        for v = 1:length(vars),
            switch vars{v}
                case 'z_correct'
                    alldata.(vars{v})(:,l) = ...
                        (dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                        dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
                case 'z_error'
                    alldata.z_error(:,l) = ...
                        (dat.(['z_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                        dat.(['z_prev' num2str(lname) 'stim__' bestmodelname]));
                case 'v_correct'
                    alldata.v_correct(:,l) = ...
                        nanmean(dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) + ...
                        dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
                case 'v_error'
                    alldata.v_error(:,l) = ...
                        (dat.(['v_prev' num2str(lname) 'resp__' bestmodelname]) - ...
                        dat.(['v_prev' num2str(lname) 'stim__' bestmodelname]));
                    
                case 'z_correct_recode'
                    alldata.z_correct_recode(:,l) = ...
                        (dat.(['z_prev' num2str(l) 'resp_correct__' bestmodelname 'recode']));
                    
                case 'z_error_recode'
                    alldata.z_error_recode(:,l) = ...
                        (dat.(['z_prev' num2str(l) 'resp_error__' bestmodelname 'recode']));
                    
                case 'v_correct_recode'
                    alldata.v_correct_recode(:,l) = ...
                        (dat.(['v_prev' num2str(l) 'resp_correct__' bestmodelname 'recode']));
                    
                case 'v_error_recode'
                    alldata.v_error_recode(:,l) = ...
                        (dat.(['v_prev' num2str(l) 'resp_error__' bestmodelname 'recode']));
                    
                
                    
            end
        end
    end
    
    % print summary
    close all;
    
    
    plotvars = {'z_correct', 'z_error', 'v_correct', 'v_error'};
    for v = 1:length(plotvars),
        subplot(2,2,v);
        plot(alldata.(plotvars{v}), alldata.([plotvars{v} '_recode']), 'o');
        grid on; r = refline(1,0); r.Color = 'k'
        title(plotvars{v});
    end
    
    print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/test_%d.pdf', d));

%     corrplot(alldata, {'z_correct', 'z_error', 'v_correct', 'v_error'}, ...
%         {'z_error_recode', 'z_correct_recode', 'v_correct_recode', 'v_error_recode'});
%   
end
end