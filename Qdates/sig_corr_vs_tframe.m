function [f,Zscores,Pvals,T_short,T_long] = sig_corr_vs_tframe(corr,tstep,days,varargin)

def_short = 14;
def_long = 180;
def_plot_title = {'RDS vs. PHQ over disjoint time-frames:',...
    'significance of correlation difference'};

p = inputParser;
addParameter(p,'short',def_short)
addParameter(p,'long',def_long)
addParameter(p,'n_subj_thresh',10)
addParameter(p,'plot_title',def_plot_title)
parse(p,varargin{:})

short_time = p.Results.short;
long_time = p.Results.long;
n_subj_thresh = p.Results.n_subj_thresh;
plot_title = p.Results.plot_title;

short_idx = tstep <= short_time;
long_idx = tstep >= long_time;

R_short = corr(short_idx);
tstep_short = tstep(short_idx);
N_short = zeros(size(R_short));
for i=1:length(tstep_short)
    N_short(i) = length(find(days <= tstep_short(i)));
end
admissable_short = N_short >= n_subj_thresh;
R_short = R_short(admissable_short);
N_short = N_short(admissable_short);

R_long = corr(long_idx);
tstep_long = tstep(long_idx);
N_long = zeros(size(R_long));
for i=1:length(tstep_long)
    N_long(i) = length(find(...
        days >= long_time & days <= tstep_long(i)));
end
admissable_long = N_long >= n_subj_thresh;
R_long = R_long(admissable_long);
N_long = N_long(admissable_long);


[RS,RL] = meshgrid(R_short,R_long);
[NS,NL] = meshgrid(N_short,N_long);

[Zscores,Pvals] = calc_rho_zscores(RS,NS,RL,NL);

[T_short,T_long] = meshgrid(tstep_short(admissable_short),...
    tstep_long(admissable_long));

f=figure;
surf(log10(T_short),log10(T_long),-log10(Pvals),RS-RL,...
    'EdgeColor','none');
title(plot_title)
xlabel('Days elapsed ($\log_{10} N$)','interpreter','latex')
ylabel('Days elapsed ($\log_{10} N$)','interpreter','latex')
zlabel('Significance ($-\log_{10}(p)$)','interpreter','latex')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'Correllation difference ($\rho_1 - \rho_2$)';
end