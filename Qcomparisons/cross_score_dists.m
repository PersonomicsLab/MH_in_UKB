function [] = cross_score_dists(varargin)
def_score_names = {'RDS','PHQ','Nr','GAD'};

p = inputParser;
addParameter(p,'Dsum',[])
addParameter(p,'see_scores',1:4)
addParameter(p,'score_names',def_score_names);
parse(p,varargin{:})

Dsum = p.Results.Dsum;
if isempty(Dsum)
    vars = load('Dsum.mat','Dscores');
    Dscores = vars.Dscores;
end
see_scores = p.Results.see_scores;
score_names = p.Results.score_names;

basic_plot_title = 'Cross-Score Distributions:';

n = length(see_scores);

for i=1:n
    for j=(i+1):n
        plot_title = {basic_plot_title,...
            [score_names{i} ' vs. ' score_names{j}]};
        
        
        val_score = Dscores(:,i);
        val_score_list = unique(val_score);
        dist_score = Dscores(:,j);
        dist_score_edges = unique(dist_score);
        
        numplts = length(val_score_list);
        dimset = split_n_dims(numplts,2);
        
        figure,
        for k=1:numplts
            scoreval = val_score_list(k);
            score1set = val_score == scoreval;
            score2set = dist_score(score1set);

            subplot(dimset(1),dimset(2),k)
            histogram(score2set,dist_score_edges)
            title(sprintf([score_names{i} '=%1d'],scoreval))
            xlabel([score_names{j} ' values'])
        end
        sgtitle(plot_title)
    end
end

end

function [dims] = split_n_dims(N,Ndims)
dims = ones(1,Ndims);
for i=0:Ndims-1
    dims(i+1) = find_first_dim(N/prod(dims),Ndims-i);
end
end

function [dim1] = find_first_dim(N,Ndims)
div_list = divisors(N);

even_split = N^(1/Ndims);

[~,nearest_fidx] = min(abs(div_list - even_split));

dim1 = div_list(nearest_fidx);
end