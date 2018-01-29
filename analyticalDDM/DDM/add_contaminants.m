% Adds uniform distribution of contaminant responses to first passage time
% densities for correct and error responses for a given set of parameters.
%
% Inputs:
%     gC_in = fptd for correct responses
%     gE_in = fptd for error responses
%     ts = times vector for fptds - dt must be fixed throughout
%     min/maxRT = min/max RT observed in the actual data & used as bounds on contaminant distribution
%     C = proportion of contaminant responses - make sure this scales with magnitude of cumulative densities

function [gC_out,gE_out] = add_contaminants(gC_in,gE_in,ts,minRT,maxRT,C)

% getting cumulative densities for later re-normalization
cum_gC=sum(gC_in);
cum_gE=sum(gE_in);

% adding uniform distrubtion of contaminants
gC_in(ts>=minRT & ts<=maxRT) = gC_in(ts>=minRT & ts<=maxRT)+((C/2)/length(find(ts>=minRT & ts<=maxRT)));
gE_in(ts>=minRT & ts<=maxRT) = gE_in(ts>=minRT & ts<=maxRT)+((C/2)/length(find(ts>=minRT & ts<=maxRT)));

% re-normalizing
if cum_gC>0   % just in case this iteration's params are so extreme as to yield no density>0 in range
    gC_out = gC_in.*(cum_gC/sum(gC_in));
else gC_out = gC_in;
end
if cum_gE>0
    gE_out = gE_in.*(cum_gE/sum(gE_in));
else gE_out = gE_in;
end

end