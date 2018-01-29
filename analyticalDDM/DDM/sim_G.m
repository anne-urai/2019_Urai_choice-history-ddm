% function for estimating probability of an incorrect response given:
% t = current time point (in seconds)
% v = drift rate
% a = boundary separation
% z = starting point
% s = noise
% Pe = probability of incorrect response

function G = sim_G(t,v,a,z,s,Pe)

tol = 10e-29; % threshold for terminating sum component of equation

sum_term=[0 0]; diff_term=[1 1];
k = 0;

while 1
    k = k+1;
    sum_term(end+1) = sum_term(end)+((2*k*sin(k*pi*z/a)*exp(-0.5*((((v.^2)/(s.^2))+((pi.^2)*(k.^2)*(s.^2)/(a.^2)))*t)))/(((v.^2)/(s.^2))+((pi.^2)*(k.^2)*(s.^2)/(a.^2))));
    diff_term(end+1) = sum_term(end)-sum_term(end-1);
    if abs(diff_term(end))<=sum_term(end-1)*tol && abs(diff_term(end-1))<=sum_term(end-2)*tol  % less than OR EQUAL TO accounts for instances where sum_term asymptotes exactly to zero - can happen in extreme cases
        break
    elseif sum_term(end-1)<0 && sum_term(end-2)<0 && abs(diff_term(end))<abs(sum_term(end-1))*tol && abs(diff_term(end-1))<abs(sum_term(end-2))*tol  % in case sum_term asymptotes to negative value - can happen in extreme cases
        break
    end
end
G = Pe-((((pi*(s.^2))/(a.^2))*exp(-(v*z/(s.^2))))*sum_term(end));