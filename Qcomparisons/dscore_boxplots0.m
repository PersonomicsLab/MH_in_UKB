function [] = dscore_boxplots(Dsum)

RecurrentMDD = logical(Dsum(:,1));
Dsum(:,1) = [];

color_labels = zeros(size(Dsum));
color_labels(RecurrentMDD,:) = 'm';
color_labels(~RecurrentMDD,:) = 'c';
color_labels = char(color_labels);

Score_groups = cell(size(Dsum));
names = {'RDS','PHQ','N','GAD'};
for i=1:length(names)
    [Score_groups{:,i}] = deal(names{i});
end

figure,
boxplot(Dsum(:),Score_groups(:),'PlotStyle','compact',...
    'ColorGroup',color_labels(:))
end