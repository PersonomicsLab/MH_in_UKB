clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')

MH = load('CCA_inputs_behavior_CCA.mat');
CCA1 = load('CCA_results_Dsum.mat');
UV = mean([CCA1.U(:,1) CCA1.V(:,1)],2);
R = corr(MH.Dep_all(:,1:4),UV);  R = R*-1;

figure
spider_plot(R',...
    'AxesLabels', {'RDS depressed mood', 'RDS disinterest', 'RDS restlessness', 'RDS tiredness'},...
    'AxesInterval', 2,...
    'AxesPrecision', 1,...
    'AxesLimits', [0.2,0.2,0.2,0.2;0.4,0.4,0.4,0.4],...
    'AxesFontSize', 8,...
    'LabelFontSize', 10);
title('Mental health CCA correlations');
print(gcf,'Fig_McDonnell','-dpng','-r300');

% Find edge numbers
load('PostCCA.mat');
good25 = [1  2  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];
good100 = [2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 45 46 48 49 50 52 53 57 58 60 63 64 93];
for n = 1:size(sig_rest_names,1)
    if isempty(cell2mat(sig_rest_names{n,2}))
        if strcmp(sig_rest_names{n,1},'PNET100')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),good100(cell2mat(sig_rest_names{n,4})),good100(cell2mat(sig_rest_names{n,5})))};
        elseif strcmp(sig_rest_names{n,1},'FNET100')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),good100(cell2mat(sig_rest_names{n,4})),good100(cell2mat(sig_rest_names{n,5})))};
        elseif strcmp(sig_rest_names{n,1},'PNET25')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),good25(cell2mat(sig_rest_names{n,4})),good100(cell2mat(sig_rest_names{n,5})))};
        elseif strcmp(sig_rest_names{n,1},'FNET25')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),good25(cell2mat(sig_rest_names{n,4})),good100(cell2mat(sig_rest_names{n,5})))};
        end
    else
        sig_rest_names{n,2} = {sprintf('%s %d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,2}))};
    end
end

% Plot resultls
R_REST(:,1) = R_REST(:,1)*-1;
R_IDP(:,1) = R_IDP(:,1) * -1;
R_REST = R_REST(sig_rest,:); R_IDP = R_IDP(sig_nonrest,:);
Rbrain = [R_REST(:,1); R_IDP(:,1)];
Nbrain = [sig_rest_names{:,2}; sig_nonrest_names];
[Rbrain,i] = sort(Rbrain); Nbrain = Nbrain(i);

figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
subtract = R_IDP(:,3); sign1 = sign(R_IDP(:,1)); sign2 = sign(R_IDP(:,3));
for n = 1:length(subtract); if sign1(n)~=sign2(n); subtract(n) = 0; end; end
barh([R_IDP(:,3) R_IDP(:,1)-subtract],'stacked');
set(gca,'ytick',1:size(R_IDP,1),'yticklabel',sig_nonrest_names(inonrest));
legend({'Correlation with UV in replication sample','Correlation with UV in exploratory CCA sample'},'Location','SouthEast')
title('A. Non-resting IDP CCA correlations');
xlabel('Cross-subject correlation with CCA score (UV)');
print(gcf,'Fig_CCA_nonrest','-dpng','-r300');

figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
I = find(R_REST(:,1)>0,1,'first');
subtract = R_REST(:,3); sign1 = sign(R_REST(:,1)); sign2 = sign(R_REST(:,3));
for n = 1:length(subtract); if sign1(n)~=sign2(n); subtract(n) = 0; end; end
subplot(1,2,1); barh(1:I-1,[R_REST(1:I-1,3) R_REST(1:I-1,1)-subtract(1:I-1)],0.6,'stacked');
set(gca,'ytick',1:I-1,'yticklabel',table2cell(sig_rest_names(irest(1:I-1),2)),'fontsize',6);
title('A. Resting IDP CCA correlations (negative)');
xlabel('Cross-subject correlation with CCA score (UV)');
subplot(1,2,2); barh(1:size(R_REST,1)-(I-1),[R_REST(I:end,3) R_REST(I:end,1)-subtract(I:end)],0.6,'stacked');
set(gca,'ytick',1:size(R_REST,1)-(I-1),'yticklabel',table2cell(sig_rest_names(irest(I:end),2)),'fontsize',6);
legend({'Correlation with UV in replication sample','Correlation with UV in exploratory CCA sample'},'Location','SouthEast')
title('B. Resting IDP CCA correlations (positive)');
xlabel('Cross-subject correlation with CCA score (UV)');
print(gcf,'Fig_CCA_rest','-dpng','-r300');
