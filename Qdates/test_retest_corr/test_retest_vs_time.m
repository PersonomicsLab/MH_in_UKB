clearvars

load('test_retest_dates_runvars.mat','Dsum12','elps_days')
labels = {'RDS','N'};

n_samp = 25000;
cmp_tst= @(x,y) corr(x,y,'type','Pearson');
method = 'perm_test';


for i=1:2
    days   = elps_days{i};
    tstep  = unique(elps_days{i});
    
    [corr_vs_days,pvals_vs_days,distributions] = ...
        time_corr_filter(Dsum12{i}(:,1),Dsum12{i}(:,2),days,...
        'n_samp',n_samp,'comparison_test',cmp_tst,'method',method,...
        'parallelize',false);
    savename = ['test_retest_vs_days_' labels{i} '.mat'];
    save(savename,'corr_vs_days','pvals_vs_days','distributions')
end

clearvars