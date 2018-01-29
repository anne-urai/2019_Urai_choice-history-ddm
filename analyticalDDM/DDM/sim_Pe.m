% function for estimating probability of an incorrect response given:
% v = drift rate
% a = boundary separation
% z = starting point
% s = noise

function Pe = sim_Pe(v,a,z,s)

Pe = (exp(-2*v*a/(s.^2))-exp(-2*v*z/(s.^2)))./(exp(-2*v*a/(s.^2))-1);

end

