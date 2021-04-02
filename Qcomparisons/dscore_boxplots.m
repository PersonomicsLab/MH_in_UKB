function [f] = dscore_boxplots(Dsum)

RecurrentMDD = logical(Dsum(:,1));
Dsum(:,1) = [];

names = {'RDS','PHQ','N','GAD'};

prob_dep_positions = 1.65:2:7.65;
no_prob_dep_positions = 2:2:8;

f=figure;
hold on
boxplot(Dsum(RecurrentMDD,:),names,'PlotStyle','compact',...
    'Colors','m','Positions',prob_dep_positions)
boxplot(Dsum(~RecurrentMDD,:),names,'PlotStyle','compact',...
    'Colors','c','Positions',no_prob_dep_positions,...
    'LabelOrientation','horizontal')
hold off
title({'Depression Score Distributions',...
    'vs. Probable MDD'})
ylabel('Score Value')

exportgraphics(f,'score_boxplots.png','Resolution',600)
end