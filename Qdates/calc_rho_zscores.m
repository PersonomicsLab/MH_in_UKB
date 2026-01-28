function [zscores,pvals] = calc_rho_zscores(r1,n1,r2,n2)
denom = sqrt(1./(n1-3) + 1./(n2-3));

zscores = (fisher_trans(r1) - fisher_trans(r2))./denom;

pvals = normcdf(zscores,'upper');
end

function [f] = fisher_trans(r)

f = 0.5*log((1+r)./(1-r));

end