function [null_dist_val,null_dist_p] = generate_null_dist(data1,data2,varargin)
rank_corr = @(x,y) corr(x,y,'type','spearman');

assert(isequal(size(data1),size(data2)),['Input variables should '...
    ' have same dimensions.'])

p = inputParser;
addParameter(p,'n_samp',1000)
addParameter(p,'comparison_test',rank_corr)
addParameter(p,'method','perm_test')
addParameter(p,'parallelize',true)
addParameter(p,'n_workers',4)
parse(p,varargin{:})

n_samp = p.Results.n_samp;
cmp_tst= p.Results.comparison_test;
method = p.Results.method;
par    = p.Results.parallelize;
n_par  = p.Results.n_workers;

if par
    if isempty(gcp('nocreate'))
        parpool(n_par)
    end
end

null_dist_val = zeros(n_samp,1);
null_dist_p = zeros(n_samp,1);

switch method
    case 'perm_test'
        for i=1:n_samp
            idx = randperm(numel(data2));
            [val,p] = cmp_tst(data1(:),data2(idx));
            null_dist_val(i) = val;
            null_dist_p(i) = p;
        end
end
end