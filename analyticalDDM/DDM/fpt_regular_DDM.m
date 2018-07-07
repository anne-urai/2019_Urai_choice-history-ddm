function [gC,gE,ts] = fpt_regular_DDM(pm,tmax)
% Peter Murphy, adapted by Anne Urai

% Get parameters
v   = pm(1);        % drift rate
Ter = pm(2);        % non-decision time
a   = pm(3);        % boundary separation
eta = pm(4);        % drift rate variability
z   = pm(5);        % starting point
s   = 0.1;          % noise

% Set timing parameters
dt = 0.01; % estimation time step (in seconds)
ts = dt:dt:(ceil((tmax-Ter)*1000)/1000);

% Estimate choice probabilities and defective cumulative probability distributions
gC=0; gE=0; % starting @ zero means ts will be equal length(diff(g))
if eta == 0 % if no drift rate variability
    Pe = sim_Pe(v,a,z,s);
    for t = dt:dt:(ceil((tmax-Ter)*1000)/1000);
        gE(end+1) = sim_G(t,v,a,z,s,Pe);
        gC(end+1) = sim_G(t,-v,a,a-z,s,1-Pe);
    end
else % if drift rate varaibility, use Gaussian quadrature to obtain G(t) by integrating over drift rate distribution
    if eta<=0.025,     n_nodes=7;
    elseif eta<=0.05,  n_nodes=8;
    elseif eta<=0.1,   n_nodes=13;   % setting number of nodes for Gaussian quadrature depending on current eta (less nodes required for lower eta)
    elseif eta<=0.175, n_nodes=15;
    elseif eta<=0.225, n_nodes=17;
    else               n_nodes=21;
    end
    [x,w] = lgwt(n_nodes,v-(4*eta),v+(4*eta));
    Fc=[]; Fe=[];
    for t = dt:dt:(ceil((tmax-Ter)*1000)/1000);
        for epsi = 1:length(x);
            Pe = sim_Pe(x(epsi),a,z,s);
            Fe(epsi,1) = sim_G(t,x(epsi),a,z,s,Pe)*(1/sqrt(2*pi*(eta.^2)))*(exp(-(((v-x(epsi)).^2)/(2*(eta.^2)))));
            Fc(epsi,1) = sim_G(t,-x(epsi),a,a-z,s,1-Pe)*(1/sqrt(2*pi*(eta.^2)))*(exp(-(((v-x(epsi)).^2)/(2*(eta.^2)))));
        end
        gE(end+1) = sum(Fe.*w);
        gC(end+1) = sum(Fc.*w);
    end
end

% Transforming cumulative distributions into normalized densities (that sum to ~1 if fully estimated)
gE = diff(gE)./(dt*1000);
gC = diff(gC)./(dt*1000);

% Getting rid of weird density values for early RTs which can arise in very extreme cases
if gE(1)>gE(2), gE(1)=gE(2); end % deals with one case in which density @ 1st time point was extremely high, and directly followed by an extremely *negative* value
if gC(1)>gC(2), gC(1)=gC(2); end % think this can happen when eta is ridiculously large and fptd is estimated for right-extreme node of drift rate distribution
gE(find(gE<0))=0; gC(find(gC<0))=0;  % replacing any negative density values (again, extremely rare) with zeros

% Time points for fpt densities plus non-decision time
ts = ts+(round(Ter*1000))/1000;

% Padding start of fpt densities with zeros
gC = [zeros(1,length([0 dt:dt:(min(ts)-dt)])) gC];
gE = [zeros(1,length([0 dt:dt:(min(ts)-dt)])) gE];
ts_extra = [0 dt:dt:(min(ts)-dt) ts];

% Interpolating fpt densities to yield ms temporal resolution within the precise time window we want (i.e. up to exactly tmax, no further)
ts = min(ts_extra):0.001:tmax;
gC = interp1(ts_extra,gC,ts,'spline');
gE = interp1(ts_extra,gE,ts,'spline');

% Setting everything before Ter and any negative values (which may result from spline interp) to zero
gC(find(ts<Ter | gC<0))=0; gE(find(ts<Ter | gE<0))=0;

