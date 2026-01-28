function [KS_results,T_results,T_stats] = dist_tests(Dsum,varargin)

p = inputParser;
addParameter(p,'display_dists',true)
parse(p,varargin{:})

display_dists = p.Results.display_dists;
names = {'RDS','PHQ','N','GAD'};

recurrent_mdd = logical(Dsum(:,1));
Dsum(:,1) = [];

KS_results = zeros(length(size(Dsum,2)),3);
T_results = zeros(length(size(Dsum,2)),2);
T_stats = cell(length(size(Dsum,2)),1);

if display_dists
    figure
end

for i=1:size(Dsum,2)
    tmp = Dsum(:,i);
    rmd_pop = tmp(recurrent_mdd);
    no_rmd_pop = tmp(~recurrent_mdd);
    [KS_results(i,1),KS_results(i,2),KS_results(i,3)] = ...
        kstest2(rmd_pop,no_rmd_pop);
    [T_results(i,1),T_results(i,2),T_stats{i}] = ...
        ttest2(rmd_pop,no_rmd_pop);
    
    if display_dists
        subplot(2,2,i)
        hold on
        histogram(rmd_pop)
        histogram(no_rmd_pop)
        hold off
        legend('Recurrent MDD','No Recurrent MDD')
        title(['Distributions of ' names{i}])
    end
end

end