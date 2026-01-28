function [] = make_figs_v2(dirname,varargin)
results_list = dir(fullfile(dirname,'*vs*.mat'));
plot_title = {'Correlation over Disjoint Times:'};

loadname = fullfile(dirname,'run_parameters.mat');
runvars = load(loadname,'num_days_elapsed','tstep');

def_position = [130,264,1010,377];
def_position2 = [125,212,573,417];

p=inputParser;
addParameter(p,'do_sig_tests',false)
addParameter(p,'plot_position',def_position)
addParameter(p,'plot2_position',def_position2)
parse(p,varargin{:})

do_sig_tests = p.Results.do_sig_tests;
position = p.Results.plot_position;
position2= p.Results.plot2_position;

for i=1:length(results_list)
    varname = results_list(i).name;
    newvars = load(fullfile(dirname,varname),'corr_vs_days');
    
    name_cell = strsplit(varname,'.');
    comp_name = name_cell{1};
    compname = strrep(comp_name,'_',' ');
    compname = strrep(compname,'vs','vs.');
    plot_title{2} = compname;
    
    f = plot_corr_vs_tframe(newvars.corr_vs_days,runvars.tstep,...
        runvars.num_days_elapsed,...
        'plot_title',plot_title);
    
    set(f, 'Position', position);
    saveloc = fullfile(dirname,[comp_name '_v2']);
    savefig(f,saveloc)
    exportgraphics(f,[saveloc '.png'],'Resolution',300)
    
    if do_sig_tests
        [g,Z,P,t_short,t_long] = sig_corr_vs_tframe(newvars.corr_vs_days,runvars.tstep,...
            runvars.num_days_elapsed,...
            'plot_title',{'Correlation Gap and Significance',...
            [plot_title{2} ' over time']});
        save(fullfile(dirname,varname),'Z','P','t_short','t_long','-append')

        set(g, 'Position', position2);
        savefig(g,[saveloc '_sig']);
        exportgraphics(g,[saveloc '_sig.png'],'Resolution',300)
    end    
end