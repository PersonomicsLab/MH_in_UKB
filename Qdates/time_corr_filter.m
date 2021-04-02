function [window_corr,window_phat,distributions] = ...
    time_corr_filter(data1,data2,time,varargin)

rank_corr = @(x,y) corr(x,y,'type','spearman');

assert(isequal(size(data1),size(data2)),['Input data1, data2, and time '...
    'should have same dimensions.'])
assert(isequal(size(data1),size(time)),['Input data1, data2, and time '...
    'should have same dimensions.'])

p = inputParser;
addParameter(p,'n_samp',1000)
addParameter(p,'comparison_test',rank_corr)
addParameter(p,'method','perm_test')
addParameter(p,'increments',unique(time))
addParameter(p,'parallelize',true)
addParameter(p,'n_workers',4)
parse(p,varargin{:})

n_samp = p.Results.n_samp;
cmp_tst= p.Results.comparison_test;
method = p.Results.method;
par    = p.Results.parallelize;
n_par  = p.Results.n_workers;
increments = p.Results.increments;

window_corr = zeros(size(increments));
window_phat = zeros(size(increments));

distributions = cell(size(increments));

for i=1:length(increments)
    time_window = (time <= increments(i));  % assumes 
    data1_recent = data1(time_window);
    data2_recent = data2(time_window);
    
    [null_dist] = generate_null_dist(data1_recent,data2_recent,...
        'n_samp',n_samp,'comparison_test',cmp_tst,'method',method,...
        'parallelize',par,'n_workers',n_par);
    distributions{i} = null_dist;
    
    window_corr(i) = cmp_tst(data1_recent,data2_recent);
    tmp_prctl = invprctile(null_dist,window_corr(i));
    tmp_phat = (100 - tmp_prctl)/100;
    
    if tmp_phat > 0
        window_phat(i) = tmp_phat;
    else
        window_phat(i) = 1e-5;
    end
end


end