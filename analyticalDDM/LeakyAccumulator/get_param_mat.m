function param_mat = get_param_mat(set_p,nsets)

param_mat=[];

if strcmp(set_p.beta,'fi')
    param_mat(1:nsets,1)=1;
elseif strcmp(set_p.beta,'fr')
    param_mat(1:nsets,1)=1:nsets;
end
if strcmp(set_p.sigmaA,'fi')
    param_mat(1:nsets,2)=1;
elseif strcmp(set_p.sigmaA,'fr')
    param_mat(1:nsets,2)=1:nsets;
end
if strcmp(set_p.theta,'fi')
    param_mat(1:nsets,3)=1;
elseif strcmp(set_p.theta,'fr')
    param_mat(1:nsets,3)=1:nsets;
end
if strcmp(set_p.leak,'no')
    param_mat(1:nsets,4)=0;
elseif strcmp(set_p.leak,'fi')
    param_mat(1:nsets,4)=1;
elseif strcmp(set_p.leak,'fr')
    param_mat(1:nsets,4)=1:nsets;
end