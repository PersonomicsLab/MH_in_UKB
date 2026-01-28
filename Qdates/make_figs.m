function [] = make_figs(dirname)
results_list = dir('*vs*.mat');
plot_title = {'Spearman Rank Corr. Coeff. over Time:'};

loadname = fullfile(dirname,'run_parameters.mat');
runvars = load(loadname,'num_days_elapsed','tstep');

position = [130,264,910,377];

for i=1:length(results_list)
    varname = results_list(i).name;
    newvars = load(varname,'corr_vs_days','distributions');
    
    name_cell = strsplit(varname,'.');
    comp_name = name_cell{1};
    compname = strrep(comp_name,'_',' ');
    compname = strrep(compname,'vs','vs.');
    plot_title{2} = compname;
    
    f = plot_corr_vs_time(newvars.corr_vs_days,runvars.tstep,...
        runvars.num_days_elapsed,...
        'dists',newvars.distributions,...
        'plot_title',plot_title,'show_null',false,...
        'plot_xlines',true,'y_label','Spearman Rank Correlation');
    
    set(f, 'Position', position);
    saveloc = fullfile(dirname,comp_name);
    savefig(f,saveloc)
    exportgraphics(f,[saveloc '.png'],'Resolution',300)
end