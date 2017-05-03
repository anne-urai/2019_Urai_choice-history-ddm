
datasets = {'RT_RDK', 'projects/0/neurodec/Data/MEG-PL', 'NatComm', 'Anke_2afc_neutral'};
datasetnames = {'RT', '2IFC',  'NatComm', 'Anke neutral', 'Anke repetitive', 'Anke alternating'};

% ============================================ %
% test correlation between dc-rep and z-rep
% ============================================ %

for d = 1:length(datasets),
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    
    % recode from stimcoding, always use 1 - 2
    dat.dc_rep__stimcodingdczprevresp = dat.dc_1__stimcodingdczprevresp - dat.dc_2__stimcodingdczprevresp;
    dat.z_rep__stimcodingdczprevresp = dat.z_1__stimcodingdczprevresp - dat.z_2__stimcodingdczprevresp;
    
    close;
    corrplot(dat, {'dc_rep__stimcodingdczprevresp', 'z_rep__stimcodingdczprevresp', ...
        'v_prevresp__regressdczprevresp', 'z_prevresp__regressdczprevresp'});
    suplabel(datasetnames{d}, 't');
    
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/biasDirectionTest_d%d.eps', d));
    
end

% ============================================ %
% is this push-pull also present at the intercept?
% ============================================ %

for d = 1:length(datasets),
    dat = readtable(sprintf('~/Data/%s/HDDM/summary/allindividualresults.csv', ...
        datasets{d}));
    
    close;
    corrplot(dat, {'v_Intercept__regressnohist', 'dc__stimcodingnohist', ...
        'z__regressnohist', 'z__stimcodingnohist'});
    suplabel(datasetnames{d}, 't');
    
    print(gcf, '-depsc', sprintf('~/Data/serialHDDM/overallBiasDirectionTest_d%d.eps', d));
    
end