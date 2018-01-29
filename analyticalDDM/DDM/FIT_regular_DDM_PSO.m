% function [pm_best,L,AIC,BIC,tr,te] = FIT_regular_DDM_PSO(subj)
%
% Fits the two-choice diffusion model to single-subject SAT data, using
% particle swarm optimzation.
% 
% Input: 'subj_in' = string variable specifying participant number
%        'r'       = string variable specifying run iteration
%        'constraints' = string variable specifying model constraints
%
% Outputs: 'pm_best' = vector of best-fitting diffusion parameters
%          'L' = likelihood value minimized by PSO procedure
%          'AIC' = Akaike Information Criterion
%          'BIC' = Bayesian Information Criterion
%          'tr' = likelihood minimum at every iteration
%          'te' = number of iterations run


function [pm_best,L,AIC,BIC,tr,te] = FIT_regular_DDM_PSO(subj,r,constraints)

% set global varaibles for passing to PSO routine
global RTc_DL RTe_DL n_miss RTc_FR RTe_FR maxRT_FR
global param_order
global C_setting

% Adding required paths
addpath(genpath('/home/pmurphy1/matlab/Imported_toolboxes'))
addpath(genpath('/home/pmurphy1/RDM_DL/Modelling/Full_fits_PSO/Full_fitting_scripts'))

% Define saving directory
loadpath = '/home/pmurphy1/RDM_DL/Modelling/Behaviour/';
savepath = '/home/pmurphy1/RDM_DL/Modelling/Full_fits_PSO/Fits/Regular_DDM/';

% Choosing parameter constraints (switch b/w 'fi' for fixed and 'fr' for 
% free) -- last two characters should be 'wC' or
% nC' to in/exclude uniform contaminant distribution
set_p.v = constraints(1:2);
set_p.Ter = constraints(4:5);
set_p.a = constraints(7:8);
set_p.eta = constraints(10:11);
set_p.C = constraints(13:14);

C_setting = set_p.C;

disp(['v=',set_p.v,', Ter=',set_p.Ter,', a=',set_p.a,', eta=',set_p.eta,', contaminants=',set_p.C])

range_p.v = [0.005 0.8];    % parameter ranges
range_p.Ter = [0.06 0.9];
range_p.a = [0.03 0.6];
range_p.eta = [0.0001 0.7];

mv_p.v = 0.05;    % maximum particle velocities
mv_p.Ter = 0.08;
mv_p.a = 0.03;
mv_p.eta = 0.03;

seeds.v.DL=[0.11 0.04]; seeds.v.FR=[0.098 0.04];  % good seed distributions for a selection of particles - [mean sd]
seeds.Ter.DL=[0.405 0.08]; seeds.Ter.FR=[0.42 0.08];
seeds.a.DL=[0.115 0.035]; seeds.a.FR=[0.145 0.035];
seeds.eta.DL=[0.16 0.035]; seeds.eta.FR=[0.1 0.035];

% Loading behavioural data
allsubj = {'06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26'};
n_blocks = 4;

subj = allsubj{subj};

% Seed random number generator
seed = round(sum(100*clock)); %never the same seed
rand('state', seed);

% RT trimming threshold for FR condition
RTstd = 4;   % first discard all trials with mean+(std*this_value) RT
RTthresh = 5;  % then discard all remaining trials with RT > this value (in secs)

% Looping through blocks & pulling out RTs/accuracy
for b = 1:n_blocks
    load([loadpath,subj,'_EEG/',subj,'_EEG_DL',num2str(b),'.mat'])
    if b==1
        ACC_DL = data_block(:,2);
        RT_DL = data_block(:,4);
    else ACC_DL(end+1:end+size(data_block,1)) = data_block(:,2);
        RT_DL(end+1:end+size(data_block,1)) = data_block(:,4);
    end
    
    load([loadpath,subj,'_EEG/',subj,'_EEG_FR',num2str(b),'.mat'])
    if b==1
        ACC_FR = data_block(:,2);
        RT_FR = data_block(:,4);
    else ACC_FR(end+1:end+size(data_block,1)) = data_block(:,2);
        RT_FR(end+1:end+size(data_block,1)) = data_block(:,4);
    end
end

% Dumping miss trials & trimming FR RT distribution
RT_DL(find(RT_DL>1.4)) = 0; % just in case any post-deadline responses happened to be registered - very rare
n_miss_DL = length(find(RT_DL==0));   % storing number of misses in DL condition, which will be included in likelihood estimates
ACC_DL = ACC_DL(find(RT_DL>0)); RT_DL = RT_DL(find(RT_DL>0));

ACC_FR = ACC_FR(find(RT_FR>0)); RT_FR = RT_FR(find(RT_FR>0)); % discard misses (there shouldn't be any in this condition though)
ACC_FR = ACC_FR(find(RT_FR<(mean(RT_FR)+(RTstd*std(RT_FR))))); RT_FR = RT_FR(find(RT_FR<(mean(RT_FR)+(RTstd*std(RT_FR))))); % discard extreme outlier trials based on STD
ACC_FR = ACC_FR(find(RT_FR<RTthresh)); RT_FR = RT_FR(find(RT_FR<RTthresh)); % discard extreme outlier trials based on absolute value

% Calculate key behaviour per condition
RTc_DL = RT_DL(find(ACC_DL==1));                 % vector of correct RTs (in seconds)
RTe_DL = RT_DL(find(ACC_DL==0));                 % vector of error RTs (in seconds)
n_miss = n_miss_DL;                           % total number of misses

RTc_FR = RT_FR(find(ACC_FR==1));                 % vector of correct RTs (in seconds)
RTe_FR = RT_FR(find(ACC_FR==0));                 % vector of error RTs (in seconds)
maxRT_FR = max([RTc_FR; RTe_FR]);        % max RT for FR condition - will only calculate first passage times up until this point


%%%%%%%%%%%%%%%%%%%%
%%% PSO settings %%%
%%%%%%%%%%%%%%%%%%%%

% Creating vector of parameter orders (0 = shared, 1 = DL only, 2 = FR only), vector of max velocities, and matrix of parameter ranges
param_order=[]; mv=[]; par_range=[];
if strcmp(set_p.v,'fi')
    param_order(end+1)=0; mv(end+1)=mv_p.v; par_range(end+1,1:2)=range_p.v;
else param_order(end+1:end+2)=[1 2]; mv(end+1:end+2)=[mv_p.v mv_p.v]; par_range(end+1:end+2,1:2)=[range_p.v; range_p.v];
end
if strcmp(set_p.Ter,'fi')
    param_order(end+1)=0; mv(end+1)=mv_p.Ter; par_range(end+1,1:2)=range_p.Ter;
else param_order(end+1:end+2)=[1 2]; mv(end+1:end+2)=[mv_p.Ter mv_p.Ter]; par_range(end+1:end+2,1:2)=[range_p.Ter; range_p.Ter];
end
if strcmp(set_p.a,'fi')
    param_order(end+1)=0; mv(end+1)=mv_p.a; par_range(end+1,1:2)=range_p.a;
else param_order(end+1:end+2)=[1 2]; mv(end+1:end+2)=[mv_p.a mv_p.a]; par_range(end+1:end+2,1:2)=[range_p.a; range_p.a];
end
if strcmp(set_p.eta,'fi')
    param_order(end+1)=0; mv(end+1)=mv_p.eta; par_range(end+1,1:2)=range_p.eta;
else param_order(end+1:end+2)=[1 2]; mv(end+1:end+2)=[mv_p.eta mv_p.eta]; par_range(end+1:end+2,1:2)=[range_p.eta; range_p.eta];
end

% Defining PSO options (see pso_Trelea_vectorized.m for details)
  P(1)=0;  P(2)=500;     P(3)=50;    P(4:13)=[1.6 1.9 0.9 0.4 400 1e-25 250 NaN 0 1];
% display  n_iterations  n_particles       acceleration, inertia, tolerance, etc

% Seeding first n particles with parameters drawn from realistic distributions
n_seeded = 25;
PSOseedValue=[];

if strcmp(set_p.v,'fi')
    PSOseedValue(1:n_seeded,end+1) = ((seeds.v.DL(1)+seeds.v.FR(1))/2)+(randn(n_seeded,1).*(sqrt(((seeds.v.DL(2).^2)+(seeds.v.FR(2).^2))/2)));
    if ~isempty(find(PSOseedValue(:,end)<range_p.v(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.v(1)),end) = range_p.v(1); end  % just in case there are any too-low values
else PSOseedValue(1:n_seeded,end+1:end+2) = [(seeds.v.DL(1)+(randn(n_seeded,1).*seeds.v.DL(2))) (seeds.v.FR(1)+(randn(n_seeded,1).*seeds.v.FR(2)))];
    if ~isempty(find(PSOseedValue(:,end-1)<range_p.v(1))), PSOseedValue(find(PSOseedValue(:,end-1)<range_p.v(1)),end-1) = range_p.v(1); end
    if ~isempty(find(PSOseedValue(:,end)<range_p.v(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.v(1)),end) = range_p.v(1); end
end
if strcmp(set_p.Ter,'fi')
    PSOseedValue(1:n_seeded,end+1) = ((seeds.Ter.DL(1)+seeds.Ter.FR(1))/2)+(randn(n_seeded,1).*(sqrt(((seeds.Ter.DL(2).^2)+(seeds.Ter.FR(2).^2))/2)));
else PSOseedValue(1:n_seeded,end+1:end+2) = [(seeds.Ter.DL(1)+(randn(n_seeded,1).*seeds.Ter.DL(2))) (seeds.Ter.FR(1)+(randn(n_seeded,1).*seeds.Ter.FR(2)))];
end
if strcmp(set_p.a,'fi')
    PSOseedValue(1:n_seeded,end+1) = ((seeds.a.DL(1)+seeds.a.FR(1))/2)+(randn(n_seeded,1).*(sqrt(((seeds.a.DL(2).^2)+(seeds.a.FR(2).^2))/2)));
    if ~isempty(find(PSOseedValue(:,end)<range_p.a(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.a(1)),end) = range_p.a(1); end  % just in case there are any too-low values
else PSOseedValue(1:n_seeded,end+1:end+2) = [(seeds.a.DL(1)+(randn(n_seeded,1).*seeds.a.DL(2))) (seeds.a.FR(1)+(randn(n_seeded,1).*seeds.a.FR(2)))];
    if ~isempty(find(PSOseedValue(:,end-1)<range_p.a(1))), PSOseedValue(find(PSOseedValue(:,end-1)<range_p.a(1)),end-1) = range_p.a(1); end
    if ~isempty(find(PSOseedValue(:,end)<range_p.a(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.a(1)),end) = range_p.a(1); end
end
if strcmp(set_p.eta,'fi')
    PSOseedValue(1:n_seeded,end+1) = ((seeds.eta.DL(1)+seeds.eta.FR(1))/2)+(randn(n_seeded,1).*(sqrt(((seeds.eta.DL(2).^2)+(seeds.eta.FR(2).^2))/2)));
    if ~isempty(find(PSOseedValue(:,end)<range_p.eta(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.eta(1)),end) = range_p.eta(1); end
elseif strcmp(set_p.eta,'fr')
    PSOseedValue(1:n_seeded,end+1:end+2) = [(seeds.eta.DL(1)+(randn(n_seeded,1).*seeds.eta.DL(2))) (seeds.eta.FR(1)+(randn(n_seeded,1).*seeds.eta.FR(2)))];
    if ~isempty(find(PSOseedValue(:,end-1)<range_p.eta(1))), PSOseedValue(find(PSOseedValue(:,end-1)<range_p.eta(1)),end-1) = range_p.eta(1); end
    if ~isempty(find(PSOseedValue(:,end)<range_p.eta(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.eta(1)),end) = range_p.eta(1); end
end

% Randomly seeding remaining particles within prespecified bounds
PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([P(3)-n_seeded,size(PSOseedValue,2)]),par_range',1);

% Running PSO routine
[output,tr,te] = pso_Trelea_vectorized('ML_from_fpts_regularDDM',length(mv),mv,par_range,0,P,'goplotpso',PSOseedValue);

% Saving output
pm_best = output(1:end-1);
L = output(end);

AIC = (2*L)+(2*length(pm_best));
BIC = (2*L)+(length(pm_best)*log(length(RTc_DL)+length(RTe_DL)+n_miss+length(RTc_FR)+length(RTe_FR)));

save([savepath,subj,'_',num2str(r),'_reg_DDM_V',set_p.v,'_TER',set_p.Ter,'_A',set_p.a,'_ETA',set_p.eta,'_',set_p.C,'.mat'],...
    'pm_best','L','AIC','BIC','tr','te','param_order','set_p','RTc_DL','RTe_DL','n_miss','RTc_FR','RTe_FR','maxRT_FR')
