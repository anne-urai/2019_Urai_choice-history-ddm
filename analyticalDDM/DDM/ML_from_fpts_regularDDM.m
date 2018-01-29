function [ML] = ML_from_fpts_regularDDM(pm)

% Retrieving global variables from initializing script
global RTc_DL RTe_DL n_miss RTc_FR RTe_FR maxRT_FR
global param_order
global C_setting

% Calculating fpt densities for each particle
for particle = 1:size(pm,1);
    
    % Pull params for each condition
    pmDL = pm(particle,find(param_order==0 | param_order==1));
    pmFR = pm(particle,find(param_order==0 | param_order==2));
    
    % Calculate first passage times per condition
    [gC_DL,gE_DL,ts_DL] = fpt_regular_DDM(pmDL,1.4);
    [gC_FR,gE_FR,ts_FR] = fpt_regular_DDM(pmFR,ceil(maxRT_FR*1000)/1000);
    
    % Creating mixture distributions of fpt densities (98%) and uniform contaminants from fastest to slowest observed RTs (2%)
    if strcmp(C_setting,'wC');
        [gC_DL,gE_DL] = add_contaminants(gC_DL,gE_DL,ts_DL,min([RTc_DL; RTe_DL]),max([RTc_DL; RTe_DL]),0.02);
        [gC_FR,gE_FR] = add_contaminants(gC_FR,gE_FR,ts_FR,min([RTc_FR; RTe_FR]),max([RTc_FR; RTe_FR]),0.02);
    end
    
    % Estimating likelihood of data given fpt densities
    sL_DL=[]; sL_FR=[];
    
    minDL = 1e-10; % values to replace zero likelihoods with (fptds here integrate to 1)
    minFR = 1e-10;
    
    if (sum(gC_DL)+sum(gE_DL))<1
        sL_DL(end+1) = n_miss*(-log(1-(sum(gC_DL)+sum(gE_DL)))); % negative log likelihood of misses
    else sL_DL(end+1) = n_miss*-log(minDL);  % in case likelihood of a miss given parameters is zero, setting to very small value
    end
    
    gC_DL(find(gC_DL<=0))=minDL; gE_DL(find(gE_DL<=0))=minDL;
    gC_FR(find(gC_FR<=0))=minFR; gE_FR(find(gE_FR<=0))=minFR;
    
    for i=1:length(RTc_DL); % negative log likelihood of correct RTs
        sL_DL(end+1) = -log(gC_DL(find(round(ts_DL.*1000)==round(RTc_DL(i)*1000))));
    end
    for i=1:length(RTe_DL); % negative log likelihood of error RTs
        sL_DL(end+1) = -log(gE_DL(find(round(ts_DL.*1000))==round(RTe_DL(i)*1000)));
    end
    
    for i=1:length(RTc_FR);
        sL_FR(end+1) = -log(gC_FR(find(round(ts_FR.*1000)==round(RTc_FR(i)*1000))));
    end
    for i=1:length(RTe_FR);
        sL_FR(end+1) = -log(gE_FR(find(round(ts_FR.*1000))==round(RTe_FR(i)*1000)));
    end
    
    ML(particle,1) = sum(sL_DL)+sum(sL_FR); % total likelihood
end


