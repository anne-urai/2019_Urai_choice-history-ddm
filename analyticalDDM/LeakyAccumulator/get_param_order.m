function param_order = get_param_order(set_p,nsets)

param_order=[];
if strcmp(set_p.beta,'fi')
    param_order(end+1)=0;
elseif strcmp(set_p.beta,'fr')
    param_order(end+1:end+size(data.fa,1))=[1:nsets];
end
if strcmp(set_p.sigmaA,'fi')
    param_order(end+1)=0;
elseif strcmp(set_p.sigmaA,'fr')
    param_order(end+1:end+size(data.fa,1))=[1:nsets];
end
if strcmp(set_p.theta,'fi')
    param_order(end+1)=0;
elseif strcmp(set_p.theta,'fr')
    param_order(end+1:end+size(data.fa,1))=[1:nsets];
end
if strcmp(set_p.leak,'fi') % this can also be set to 'no', in which case leak will be fixed to 0
    param_order(end+1)=0;
elseif strcmp(set_p.leak,'fr')
    param_order(end+1:end+size(data.fa,1))=[1:nsets];
end