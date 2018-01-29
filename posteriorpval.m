function p = posteriorpval(dat1, dat2)

% Code to fit the history-dependent drift diffusion models described in
% Urai AE, Gee JW de, Donner TH (2018) Choice history biases subsequent evidence accumulation. bioRxiv:251595
%
% MIT License
% Copyright (c) Anne Urai, 2018
% anne.urai@gmail.com

overlap = mean((dat1 - dat2) > 0);
p = min([overlap, 1-overlap]);

end
