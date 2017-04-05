function y = logistic(p, x)
% Parameters: p(1) bias
%             p(2) slopw
%             p(3) lapse rate
%             x   intensity values.

% include a lapse rate, see Wichmann and Hill parameterisation
y =  p(3)+(1-p(3)-p(3)) * (1./(1+exp(- (p(2).*x +p(1) ))));

end