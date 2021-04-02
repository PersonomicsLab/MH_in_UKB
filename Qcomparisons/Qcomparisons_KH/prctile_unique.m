function [unq_ptl,cutoffs] = prctile_unique(X)

[vals] = prctile(X,[1:100]);
[cutoffs,unq_ptl] = unique(vals);

end


