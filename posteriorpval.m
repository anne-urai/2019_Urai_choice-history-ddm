function p = posteriorpval(dat1, dat2)

overlap = mean((dat1 - dat2) > 0);
p = min([overlap, 1-overlap]);

end
