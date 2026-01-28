function [zstat pval df] = steiger_z(r12,r13,r23, N)
% This function computes the z test statistic and two-tailed p-value for 
% a comparison between two dependent correlations involving variables
% X1,X2,and Y using Steiger's z-test  
% inputs: r12 - corr(X1,Y), r13 - corr(X2,Y), r23 - corr(X1,X2), N - sample
% size. 
% outputs: zstat - z statistic, pval - two-tailed p-value, df - degrees of
% freedom. 
% by Joseph Griffis, 2018, Washington University in St. Louis

% Fisher z transform
z12 = fisherz(r12);
z13 = fisherz(r13);

% Compute z-statistic
zrbar = 0.5*(z12+z13);
rbar = fisherz_to_r(zrbar);
rbar2 = rbar^2;
sbar = [(r23)*(1-rbar2-rbar2)-0.5*(rbar2)*(1-rbar2-rbar2-r23^2)]./[(1-rbar2)*(1-rbar2)];
df = N-3;
zstat = [sqrt(df)*(z12-z13)]./[sqrt(2-2*sbar)];

% get p-value
pval = 1-cdf('normal', abs(zstat), 0 ,1);
pval = pval.*2;

end