function [zstat pval df] = steiger_z2(x1,x2,y, type)
% This function computes the z test statistic and two-tailed p-value for 
% a comparison between two dependent correlations involving variables
% X1,X2,and Y using Steiger's z-test  
% inputs: r12 - corr(X1,Y), r13 - corr(X2,Y), r23 - corr(X1,X2), N - sample
% size. 
% outputs: zstat - z statistic, pval - two-tailed p-value, df - degrees of
% freedom. 
% by Joseph Griffis, 2018, Washington University in St. Louis

if strcmp(type,'rank')==1
    r12 = corr(x1,y, 'type', 'spearman');
    r13 = corr(x2,y, 'type', 'spearman');
    r23 = corr(x1,x2, 'type', 'spearman');
elseif strcmp(type, 'linear')==1
    r12 = corr(x1,y);
    r13 = corr(x2,y);
    r23 = corr(x1,x2);
end

N = numel(find(isnan(y)==0));

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