preamble = fullfile(filesep,'scratch','tyo8teasley','Qdates');
% preamble = fullfile(filesep,'scratch','janine','MentalHealthInUKB','Qdates');
% preamble = fullfile('');
addpath(preamble)
loadname = fullfile(preamble,'Dates_and_scores.mat');
vars = load(loadname,'Dsum','elps_days','summary_labels');

metric_names = vars.summary_labels(2:end);
vars.Dsum(:,1) = [];
n_metrics = 4;

assert(size(vars.Dsum,2)==n_metrics,'Not all metrics are present.')
clear n_metrics;

Dsum = vars.Dsum;
num_days_elapsed = abs(vars.elps_days);
clear vars

n_samp = 25000;
cmp_tst= @(x,y) corr(x,y,'type','spearman');
method = 'perm_test';
par    = false;
tstep  = unique(num_days_elapsed);

savedir = ['metrics_thru_time_' datestr(datetime,'mmmdd')];
saveloc = fullfile(preamble,savedir);
if ~exist(saveloc,'dir')
    mkdir(saveloc)
end

cd(saveloc)

savename = fullfile(saveloc,'run_parameters.mat');
save(savename,'Dsum','num_days_elapsed','n_samp',...
    'cmp_tst','method','par','tstep')

for i=1:size(Dsum,2)
    for j=i:size(Dsum,2)
        if j~=i
            [corr_vs_days,pvals_vs_days,distributions] = ...
                time_corr_filter(Dsum(:,i),Dsum(:,j),...
                num_days_elapsed,'n_samp',n_samp,...
                'comparison_test',cmp_tst,'method',method,...
                'parallelize',par,...
                'increments',tstep);
            savename_ij = [metric_names{i} '_vs_' metric_names{j} '.mat'];
            saveloc_ij = fullfile(saveloc,savename_ij);
            save(saveloc_ij,'corr_vs_days','pvals_vs_days','distributions')
        end
    end
end