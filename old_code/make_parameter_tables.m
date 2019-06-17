function make_parameter_tabless

addpath(genpath('~/code/Tools'));
warning off; close all;
global datasets datasetnames mypath

mdls = {'z_prevresp', 'dc_prevresp', ...
    'dc_z_prevresp', 'nohist'};
tmpNames = {'deGee2014', 'deGee2017', 'Murphy2014', 'Braun', 'Urai2017', 'Urai'};

for m = 1:length(mdls),
    clear p allparams;
    
    for d = 1:length(datasets),
        
        if exist(sprintf('%s/%s/stimcoding_%s/results-combined.csv',  mypath, datasets{d}, mdls{m}), 'file'),
            p = readtable(sprintf('%s/%s/stimcoding_%s/results-combined.csv',  mypath, datasets{d}, mdls{m}));
            
            % if there are several levels of drift rate (depending on
            % stimulus difficulty), display only one in the table and point
            % to the sanity check supplementary figure
            if any(ismember(p.Var1, 'v(0.09)')),
                p.Var1(ismember(p.Var1, 'v(0.09)')) = {'v'};
            elseif any(ismember(p.Var1, 'v(0.1)')),
                p.Var1(ismember(p.Var1, 'v(0.1)')) = {'v'};
            end
            
            % ONLY GRAB THE GROUP-LEVEL PARAMETERS
            pickTheseParams = {'a', 'v', 't', 'sv', 'sz', 'z(1)', 'z(-1)', 'dc(1)', 'dc(-1)', 'dc', 'z'};
            params{d} = p(ismember(p.Var1, pickTheseParams), 1:2);
            
            % rename some things for easier concatenation
            params{d}.Properties.VariableNames{'mean'} = tmpNames{d};
            params{d}.Properties.RowNames = params{d}.Var1;
            params{d}.Var1 = [];
        end
    end
    
    %% concatenate, write to a nice-looking excel table?
    
    allparams = cat(2, params{~cellfun(@isempty, params)});
    writetable(allparams, sprintf('~/Data/serialHDDM/params_%s.xls', mdls{m}), 'WriteRowNames',true, 'WriteVariableNames', true);

    
end
end