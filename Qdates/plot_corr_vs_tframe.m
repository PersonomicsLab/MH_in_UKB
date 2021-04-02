function [f] = plot_corr_vs_tframe(corr,tstep,days,varargin)
def_short = 14;
def_long = 93;
def_plot_title = 'RDS vs. PHQ over disjoint time-frames';

p = inputParser;
addParameter(p,'short',def_short)
addParameter(p,'long',def_long)
addParameter(p,'do_middle',true)
addParameter(p,'tscale','log10lin')
addParameter(p,'plot_title',def_plot_title)
parse(p,varargin{:})

short_time = p.Results.short;
long_time = p.Results.long;
do_middle = p.Results.do_middle;
tscale = p.Results.tscale;
plot_title = p.Results.plot_title;

short_idx = tstep <= short_time;
long_idx = tstep >= long_time;
middle_idx = ~short_idx & ~long_idx;

leg_entries = {'"Short" timeframe','','"Long" timeframe',''};

plot_corr_vs_time(corr(short_idx),tstep(short_idx),...
    days,'use_dists',false,'tscale',tscale,...
    'linevals',short_time);
f = plot_corr_vs_time(corr(long_idx),tstep(long_idx),...
    days(days >= long_time),'tscale',tscale,...
    'use_dists',false,'new_fig',false,'linevals',long_time,...
    'plot_title',plot_title);

if do_middle
    plot_corr_vs_time(corr(middle_idx),tstep(middle_idx),...
        days(days>=short_time & days<=long_time),'plot_xlines',false,...
        'use_dists',false,'new_fig',false,'tscale',tscale,...
        'plot_title',plot_title);
    
    leg_entries{5} = '"Middle" time-frame';
end

legend(leg_entries,'location','southwest')
end