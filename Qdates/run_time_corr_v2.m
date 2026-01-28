% preamble = fullfile(filesep,'scratch','tyo8teasley','Qdates');
% preamble = fullfile(filesep,'scratch','janine','MentalHealthInUKB','Qdates');
preamble = fullfile('');
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

cmp_tst= @(x,y) corr(x,y,'type','spearman');
tstep  = unique(num_days_elapsed);

savedir = ['metrics_thru_time_' datestr(datetime,'mmmdd')];
saveloc = fullfile(preamble,savedir);
if ~exist(saveloc,'dir')
    mkdir(saveloc)
end

savename = fullfile(saveloc,'run_parameters.mat');
save(savename,'Dsum','num_days_elapsed','cmp_tst','tstep')

for i=1:size(Dsum,2)
    for j=i:size(Dsum,2)
        if j~=i
            [corr_vals,samp_size,tstep_new] = ...
                time_corr_filter_v2(Dsum(:,i),Dsum(:,j),...
                num_days_elapsed,'increments',tstep,...
                'comparison_test',cmp_tst);
            savename_ij = [metric_names{i} '_vs_' metric_names{j} '.mat'];
            saveloc_ij = fullfile(saveloc,savename_ij);
            save(saveloc_ij,'corr_vals','samp_size','tstep_new')
        end
    end
end