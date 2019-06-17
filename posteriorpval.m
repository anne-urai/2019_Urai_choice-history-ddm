function p = posteriorpval(dat1, dat2)

% Code to fit the history-dependent drift diffusion models as described in
% Urai AE, de Gee JW, Tsetsos K, Donner TH (2019) Choice history biases subsequent evidence accumulation. eLife, in press.
%
% MIT License
% Copyright (c) Anne Urai, 2019
% anne.urai@gmail.com

overlap = mean((dat1 - dat2) > 0);
p = min([overlap, 1-overlap]);

end
