function [f] = plot_corr_vs_time(corr,tstep,days,varargin)

assert(isequal(numel(corr),numel(tstep)),...
    'Correlation and time-step variables have different numbers of elements')

p=inputParser;
addParameter(p,'use_dists',true)
addParameter(p,'dists',zeros(size(corr)))

addParameter(p,'n_subj_thresh',10)
addParameter(p,'plot_title','Placeholder Title')
addParameter(p,'new_fig',true)
addParameter(p,'show_null',false)
addParameter(p,'plot_xlines',true)
addParameter(p,'linevals',[10 25 100])
addParameter(p,'y_label','Correlation')
addParameter(p,'tscale','loglin')
parse(p,varargin{:})

use_dists = p.Results.use_dists;
dists = p.Results.use_dists;

n_subj_thresh = p.Results.n_subj_thresh;
plot_title = p.Results.plot_title;
new_fig = p.Results.new_fig;
show_null = p.Results.show_null;
plot_xlines = p.Results.plot_xlines;
linevals = p.Results.linevals;
y_label = p.Results.y_label;
tscale = p.Results.tscale;

n_subj_per_corr = zeros(size(corr));
std_err = zeros(size(corr));
mean_rand = zeros(size(corr));
for i=1:length(tstep)
    n_subj_per_corr(i) = length(find(days <= tstep(i)));
    if use_dists
        std_err(i) = sqrt(var(dists{i}));
        mean_rand(i) = mean(dists{i});
    else
        denom = max(eps,sqrt(n_subj_per_corr(i)-1));
        std_err(i) = 0.6325/denom;
        mean_rand(i) = 0;
    end
end

admissable = (n_subj_per_corr >= n_subj_thresh);

corr_new = corr(admissable);
tstep_new = tstep(admissable);
corr_std_err = std_err(admissable);

mstyle = '-d';
mcol = 'm';

switch tscale
    case 'lin'
        xvals = tstep_new;
        x_label = 'Days Elapsed';
    case 'log10lin'
        xvals = log10(tstep_new);
        x_label = 'Days Elapsed ($\log_{10} N$)';
    case 'loglin'
        xvals = log(tstep_new);
        x_label = 'Days Elapsed ($\log N$)';
end

if new_fig
    f=figure;
else
    f=gcf;
end

hold on
if show_null
    mean_rand_new = mean_rand(admissable);
    errorbar(xvals,corr_new,corr_std_err,mstyle,...
        'MarkerFaceColor',mcol)
    errorbar(xvals,mean_rand_new,corr_std_err,mstyle,...
        'MarkerFaceColor',mcol)
    legend('Measured','Null')
else
    errorbar(xvals,corr_new,corr_std_err,mstyle,...
        'MarkerFaceColor',mcol)
end
if plot_xlines
    plot_nearest_xlines(xvals,tstep_new,linevals);
end
hold off
title(plot_title)
xlabel(x_label,'interpreter','latex')
ylabel(y_label,'interpreter','latex')

end

function [] = plot_nearest_xlines(xvals,tstep,days_array)
for i=1:length(days_array)
    plot_nearest_xline(xvals,tstep,days_array(i))
end
end

function [] = plot_nearest_xline(xvals,tstep,n_days)
[~,idx] = min(abs(tstep - n_days));

lineval = xvals(idx);
label = sprintf('%1d days elapsed',tstep(idx));

xline(lineval,'r--',label)
end