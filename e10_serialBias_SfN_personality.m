function e10_serialBias_SfN_personality
% from the csv table, make an overview of repetition behaviour

% get a huge list with values for each participant
% can then work with this dataframe

clear; close all; clc;
addpath(genpath('~/code/Tools'));

global mypath
% neutral vs biased plots
datasets = {'MEG', 'MEG_MEGsessions'};
datasetnames = {{'2IFC, replication'}, ...
    {'2IFC, replication', 'MEG sessions'}};

% ============================================ %
% ONE LARGE PLOT WITH PANEL FOR EACH DATASET
% ============================================ %

close all;
for d = length(datasets):-1:1
    disp(datasets{d});
    
    results = readtable(sprintf('%s/summary/%s/allindividualresults.csv', mypath, datasets{d}));
    results = results(results.session == 0, :);
    
    allresults = {};
    switch datasets{d}
        case {'Anke_2afc_sequential'}
            
            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY
            
            % use the stimcoding difference only from alternating
            allresults{1}.z_prevresp        = results.z_1_neu__stimcodingdczprevresp - results.z_2_neu__stimcodingdczprevresp;
            allresults{1}.v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_neu__stimcodingdczprevresp;
            allresults{1}.criterionshift    = results.repetition_neutral;
            allresults{1}.subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');
            
            allresults{2}.z_prevresp        = results.z_1_alt__stimcodingdczprevresp - results.z_2_alt__stimcodingdczprevresp;
            allresults{2}.v_prevresp        = results.dc_1_alt__stimcodingdczprevresp - results.dc_2_alt__stimcodingdczprevresp;
            allresults{2}.criterionshift    = results.repetition_alternating;
            allresults{2}.subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');
            
            allresults{3}.z_prevresp        = results.z_1_rep__stimcodingdczprevresp - results.z_2_rep__stimcodingdczprevresp;
            allresults{3}.v_prevresp        = results.dc_1_neu__stimcodingdczprevresp - results.dc_2_rep__stimcodingdczprevresp;
            allresults{3}.criterionshift    = results.repetition_repetitive;
            allresults{3}.subjnr            = results.subjnr;
            alltitles{3}                    = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');
            
        case {'Bharath_fMRI', 'Anke_MEG', 'Anke_merged'}
            
            % DOES THIS DATASET HAVE MULTIPLE TRANSITION PROBABILITIES?
            % THEN PLOT THESE SEPARATELY
            
            % use the stimcoding difference only from alternating
            allresults{1}.z_prevresp        = results.z_1_0_50_0__stimcodingdczprevresp - results.z_2_0_50_0__stimcodingdczprevresp;
            allresults{1}.v_prevresp        = results.dc_1_0_50_0__stimcodingdczprevresp - results.dc_2_0_50_0__stimcodingdczprevresp;
            allresults{1}.criterionshift    = results.repetition_neutral;
            allresults{1}.subjnr            = results.subjnr;
            alltitles{1}                    = cat(2, datasetnames{d}{1}, ' - ', 'Neutral');
            
            allresults{2}.z_prevresp        = results.z_1_0_20_0__stimcodingdczprevresp - results.z_2_0_20_0__stimcodingdczprevresp;
            allresults{2}.v_prevresp        = results.dc_1_0_20_0__stimcodingdczprevresp - results.dc_2_0_20_0__stimcodingdczprevresp;
            allresults{2}.criterionshift    = results.repetition_alternating;
            allresults{2}.subjnr            = results.subjnr;
            alltitles{2}                    = cat(2, datasetnames{d}{1}, ' - ', 'Alternating');
            
            allresults{3}.z_prevresp        = results.z_1_0_80_0__stimcodingdczprevresp - results.z_2_0_80_0__stimcodingdczprevresp;
            allresults{3}.v_prevresp        = results.dc_1_0_80_0__stimcodingdczprevresp - results.dc_2_0_80_0__stimcodingdczprevresp;
            allresults{3}.criterionshift    = results.repetition_repetitive;
            allresults{3}.subjnr            = results.subjnr;
            alltitles{3}                    = cat(2, datasetnames{d}{1}, ' - ', 'Repetitive');
            
        otherwise
            
            try
                % use the stimcoding difference
                results.z_prevresp = ...
                    results.z_1__stimcodingdczprevresp - results.z_2__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1__stimcodingdczprevresp - results.dc_2__stimcodingdczprevresp;
            catch
                results.z_prevresp = ...
                    results.z_1_0__stimcodingdczprevresp - results.z_2_0__stimcodingdczprevresp;
                results.v_prevresp = ...
                    results.dc_1_0__stimcodingdczprevresp - results.dc_2_0__stimcodingdczprevresp;
            end
            
            results.criterionshift = results.repetition;
            allresults{1} = results;
            
            alltitles{1} = datasetnames{d}{1}; % use only the dataset title
    end
    
    
    for resultsSet = 1:length(allresults),
        
        close all;
        corrplot(allresults{resultsSet}, {'repetition', 'z_prevresp', 'v_prevresp'}, ...
            {'BIS', 'BAS', 'PSWQ', 'AQ'});
        
        print(gcf, '-dpdf', sprintf('~/Data/serialHDDM/personality_d%d_s%d.pdf', d, resultsSet));
    end
end

end
