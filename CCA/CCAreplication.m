clear all; close all; clc

addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/FSLNets')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/FSL')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/')
INPUT = 'Subjects_CCAnew.csv';
INPUT2 = 'Subjects_effectsize.csv';

% Load data
S1 = load(sprintf('%s/%s','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/SubjectSplits/',INPUT));
S2 = load(sprintf('%s/%s','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/SubjectSplits/',INPUT2));
load('Eigs_CCAnew.mat');

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s.mat','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/',INPUT2(10:end-4)));
conf2 = nets_demean(conf);
Pconf2 = pinv(conf2); clear conf


%% Resting state IDPs
fprintf('Processing resting state IDPs\n');
load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/IDPs.mat');
IDP2 = IDP; clear IDP;
%IDP2 = load('/Users/janinebijsterbosch/Box/00_WashU/Data/UKB/IDP/IDPs_rest_40k.txt');
%subs = IDP2(:,1); IDP2 = IDP2(:,2:end);
[~,s,~] = intersect(subs,S2); IDP2 = IDP2(s,:);
IDP2 = nets_demean(IDP2);
IDP2 = nets_demean(IDP2-conf2*(Pconf2*IDP2));
SCORE_PROJECT = IDP2*REST_COEFF;
REST_OUT2 = SCORE_PROJECT(:,1:Nkeep_rest);  clear SCORE_PROJECT

%% Load data
DATA = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/IDP_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S2); DATA2 = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/ExtractVariables/vars.mat','IDP_nonrest');

%% "impute" missing  data - actually this avoids any imputation
fprintf('Processing IDPs\n');
varsd = table2array(DATA2);
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
varsd(:,ID) = [];
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf2(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd = knnimpute(varsd);

%% task
fprintf('Dimensionality reduction task\n');
task1 = strncmp('Median BOLD',IDP_nonrest(:,2),11); task1 = find(task1==1);
task2 = strncmp('Median z-stat',IDP_nonrest(:,2),13); task2 = find(task2==1);
task3 = strncmp('90th percentile',IDP_nonrest(:,2),15); task3 = find(task3==1);
task = sort([task1; task2; task3]); clear task1 task2 task3
varsd2_task = varsd(:,task);
SCORE_PROJECT = varsd2_task * TASK_COEFF;
TASK_OUT2 = SCORE_PROJECT(:,1:Nkeep_task); clear SCORE_PROJECT

%% Dimensionality reduction on structural IDPs
fprintf('Dimensionality reduction structural\n');
nottask = setdiff(1:size(varsd,2),task);
varsd2_nottask = varsd(:,nottask);
SCORE_PROJECT = varsd2_nottask * STRUCT_COEFF;
STRUCT_OUT2 = SCORE_PROJECT(:,1:Nkeep_struct); clear SCORE_PROJECT

%% Multiply by CCA weights
CCA = load('CCA_results_CCAnew.mat');
MH2 = load('CCArepliation_behavior_effectsize.mat');
Vrepl = [REST_OUT2 STRUCT_OUT2 TASK_OUT2] * CCA.A;
Urepl = MH2.Dsum * CCA.B;
[r,p] = corr(Vrepl,Urepl)
p = p*5

%% Average UV for follow-up tests
UV = mean([CCA.U(:,1) CCA.V(:,1)],2);
UV2 = mean([Urepl(:,1) Vrepl(:,1)],2);

%% Correlation with eigenvectors rather than IDPs (for reviewer)
TEST = zeros(126,2);
load('CCA_inputs_brain_CCAnew.mat');
TEST(:,1) = corr(UV,[REST_OUT STRUCT_OUT TASK_OUT])';
TEST(:,2) = corr(UV2,[REST_OUT2 STRUCT_OUT2 TASK_OUT2])';
figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
barh(TEST(:,1),0.5,'FaceColor','r')
hold on
barh(TEST(:,2),0.25,'FaceColor','b')
legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
title('Eigenvector CCA correlations');
xlabel('Cross-subject correlation with CCA score (UV)'); ylabel('Eigenvectors');
print(gcf,'Fig_for_reviewer_eigenvector_UV','-dpng','-r300');

%% Behavioral correlations
MH = load('CCA_inputs_behavior_CCAnew.mat');
R = corr(MH.Dsum,UV);
R2 = corr(MH2.Dsum,UV2);
figure
spider_plot([R2 R]',...
    'AxesLabels', {'Depression status','RDS-4','PHQ-9','N-12','GAD-7'},...
    'AxesInterval', 4,...
    'AxesPrecision', 2,...
    'AxesLimits', [-1,-1,-1,-1,-1;0,0,0,0,0],...
    'Color', [0, 0, 1; 1, 0, 0],...
    'AxesFontSize', 8,...
    'LabelFontSize', 10);
legend({'Correlation with UV in replication sample','Correlation with UV in exploratory CCA sample'},'Location','SouthEast')
title('B. Mental health CCA correlations');
print(gcf,'Fig_CCA_mh','-dpng','-r300');


%% IDP correlations
load('PostCCA.mat');

r = corr(varsd,UV2(:,1),'rows','pairwise');
r = r(Iidp); R_IDP(:,3) = r; clear r
r = corr(IDP2,UV2(:,1)); R_REST(:,3) = r; clear r

% Find edge numbers
good25 = [1  2  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];
good100 = [2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 45 46 48 49 50 52 53 57 58 60 63 64 93];
I = find(strcmp(sig_rest_names.IDP_type,'AMP25'));
I2 = find(strcmp(sig_rest_names.IDP_type,'FNET25'));
I3 = find(strcmp(sig_rest_names.IDP_type,'PNET25'));
I = sort([I;I2;I3]); clear I2 I3
sig_rest_names = sig_rest_names(I,:);
sig_rest = sig_rest(I,:);
for n = 1:size(sig_rest_names,1)
        if strcmp(sig_rest_names{n,1},'AMP25')
            sig_rest_names{n,2} = {sprintf('%s %d',cell2mat(sig_rest_names{n,1}),find(good25==cell2mat(sig_rest_names{n,2})))};
        elseif strcmp(sig_rest_names{n,1},'AMP100')
            sig_rest_names{n,2} = {sprintf('%s %d',cell2mat(sig_rest_names{n,1}),find(good100==cell2mat(sig_rest_names{n,2})))};
        elseif strcmp(sig_rest_names{n,1},'PNET100')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,4}),cell2mat(sig_rest_names{n,5}))};
        elseif strcmp(sig_rest_names{n,1},'FNET100')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,4}),cell2mat(sig_rest_names{n,5}))};
        elseif strcmp(sig_rest_names{n,1},'PNET25')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,4}),cell2mat(sig_rest_names{n,5}))};
        elseif strcmp(sig_rest_names{n,1},'FNET25')
            sig_rest_names{n,2} = {sprintf('%s %d-%d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,4}),cell2mat(sig_rest_names{n,5}))};
        end
%     else
%         sig_rest_names{n,2} = {sprintf('%s %d',cell2mat(sig_rest_names{n,1}),cell2mat(sig_rest_names{n,2}))};
%     end
end

% Plot resultls
R_REST(:,1) = R_REST(:,1);
R_IDP(:,1) = R_IDP(:,1);
R_REST = R_REST(sig_rest,:); [~,irest] = sort(R_REST(:,1)); R_REST = R_REST(irest,:);
R_IDP = R_IDP(sig_nonrest,:); [~,inonrest] = sort(R_IDP(:,1)); R_IDP = R_IDP(inonrest,:);

figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
barh(R_IDP(:,1),0.5,'FaceColor','r')
hold on
barh(R_IDP(:,3),0.25,'FaceColor','b')
set(gca,'ytick',1:size(R_IDP,1),'yticklabel',sig_nonrest_names(inonrest));
legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
title('A. Non-resting IDP CCA correlations');
xlabel('Cross-subject correlation with CCA score (UV)');
%print(gcf,'Fig_CCA_nonrest','-dpng','-r300');

figure; set(gcf,'Position',[100 100 1300 950],'PaperPositionMode','auto')
I = find(R_REST(:,1)>0,1,'first');
subplot(1,2,1); 
barh(1:I-1,R_REST(1:I-1,1) ,0.5,'FaceColor','r')
hold on
barh(1:I-1,R_REST(1:I-1,3) ,0.25,'FaceColor','b')
set(gca,'ytick',1:I-1,'yticklabel',table2cell(sig_rest_names(irest(1:I-1),2)),'fontsize',6);
title('A. Resting IDP CCA correlations (negative)');
xlabel('Cross-subject correlation with CCA score (UV)');
subplot(1,2,2); 
barh(1:size(R_REST,1)-(I-1),R_REST(I:end,1) ,0.5,'FaceColor','r')
hold on
barh(1:size(R_REST,1)-(I-1),R_REST(I:end,3) ,0.25,'FaceColor','b')
set(gca,'ytick',1:size(R_REST,1)-(I-1),'yticklabel',table2cell(sig_rest_names(irest(I:end),2)),'fontsize',6);
legend({'Correlation with UV in exploratory CCA sample','Correlation with UV in replication sample'},'Location','SouthEast')
title('B. Resting IDP CCA correlations (positive)');
xlabel('Cross-subject correlation with CCA score (UV)');
print(gcf,'Fig_CCA_rest','-dpng','-r300');

