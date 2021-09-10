clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSL')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'Subjects_CCA.csv';
INPUT2 = 'Subjects_effectsize.csv';

% Load data
S1 = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
S2 = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT2));
variance_threshold = 0.5; % keep number of PCs that explain at least 50% of variance
NkeepMax = floor(length(S1)/10/2);

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s.mat','/scratch/janine/MentalHealthInUKB/Data/',INPUT(10:end-4)));
conf1 = nets_demean(conf);
Pconf1 = pinv(conf1); clear conf
load(sprintf('%s/Confounds_Subjects_%s.mat','/scratch/janine/MentalHealthInUKB/Data/',INPUT2(10:end-4)));
conf2 = nets_demean(conf);
Pconf2 = pinv(conf2); clear conf

%% Resting state IDPs
fprintf('Processing resting state IDPs\n');
load('/scratch/janine/MentalHealthInUKB/Data/IDPs.mat');
[~,s,~] = intersect(subs,S1); IDP1 = IDP(s,:)';
IDP1 = nets_demean(IDP1,2);
IDP1 = nets_demean((IDP1'-conf1*(Pconf1*IDP1'))',2);
C = IDP1' * IDP1;
[V,D] = eig(C); Vrest = V; Drest = D;
[~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
REST_OUT = V(:,1:Nkeep);

[~,s,~] = intersect(subs,S2); IDP2 = IDP(s,:)';
IDP2 = nets_demean(IDP2,2);
IDP2 = nets_demean((IDP2'-conf2*(Pconf2*IDP2'))',2);
V = calc_eig_proj(IDP2',IDP1',Vrest,Drest);
V = V(:,size(V,2):-1:1);
REST2_OUT = V(:,1:Nkeep);

%% Load data
DATA = readtable('/scratch/janine/MentalHealthInUKB/Data/IDP_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S1); DATA1 = DATA(s,:);
[~,s,~] = intersect(table2array(DATA(:,1)),S2); DATA2 = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/scratch/janine/MentalHealthInUKB/Data/ExtractVariables/vars.mat','IDP_nonrest');

%% "impute" missing  data - actually this avoids any imputation
fprintf('Processing IDPs\n');
varsd = table2array(DATA1);
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
varsd(:,ID) = [];
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf1(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsdCOV = zeros(size(varsd,1));
for i = 1:size(varsd,1) % estimate "pairwise" covariance, ignoring missing data
    for j = 1:size(varsd,1)
        grot = varsd([i j],:); 
        grot = cov(grot(:,sum(isnan(grot))==0)'); 
        varsdCOV(i,j) = grot(1,2);
    end
end
varsdCOV2 = nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[V,D] = eig(varsdCOV2); Vnonrest = V; Dnonrest = D;
[~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
IDP_OUT = V(:,1:Nkeep);

varsd1 = knnimpute(varsd);
varsd1 = varsd1';
C = varsd1' * varsd1;
[V1,D1] = eig(C); Vnonrest1 = V1; Dnonrest1 = D1;
[~,inds] = sort(diag(D1), 'descend'); D1 = D1(inds,inds); V1 = V1(:,inds);
Vexp = diag(D1)/sum(diag(D1)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
IDP_OUT_NEW = V1(:,1:Nkeep);

varsd = table2array(DATA2);
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
varsd(:,ID) = [];
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf1(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd = knnimpute(varsd);
V = calc_eig_proj(varsd,varsd1',Vnonrest1,Dnonrest1);
V = V(:,size(V,2):-1:1);
IDP2_OUT = V(:,1:Nkeep);

%% Behavioral data
% Load data
S = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
MH = readtable('/scratch/janine/MentalHealthInUKB/Data/MH_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(MH(:,1)),S); MH = MH(s,:);
H = get_UKB_headers(MH);
MH = standardizeMissing(MH,-3); MH = standardizeMissing(MH,-1); MH = standardizeMissing(MH,-818);

% Full data 
idx = ismissing(MH(:,{'x4609_2_0','x5375_2_0'}));
MH{:,{'x4609_2_0','x5375_2_0'}}(idx) = 0;
Dall = table2array(MH)';
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
Dall(ID,:) = [];
Dall = nets_demean(Dall,2);
Dall = nets_demean((Dall'-conf2*(Pconf2*Dall'))',2);
Dall = Dall';

% Recurrrent MDD
n1 = strfind(H,'4598-2.0'); n1 = find(~cellfun(@isempty,n1));
EverDep = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'4631-2.0'); n1 = find(~cellfun(@isempty,n1));
EverUnenth = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'4609-2.0'); n1 = find(~cellfun(@isempty,n1));
DurDep = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'5375-2.0'); n1 = find(~cellfun(@isempty,n1));
DurUnenth = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2090-2.0'); n1 = find(~cellfun(@isempty,n1));
SeenGP = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2100-2.0'); n1 = find(~cellfun(@isempty,n1));
SeenPsych = table2array(MH(:,n1)); clear n1
RecurrentMDD = zeros(size(S));
for n = 1:length(S)
    if EverDep(n) == 1 || EverUnenth(n) == 1
        if DurDep(n)>1 || DurUnenth(n)>1
            if SeenGP(n) == 1 || SeenPsych(n) == 1
                RecurrentMDD(n) = 1;
            end
        end
    end
end

% Recent MDD symptoms
n1 = strfind(H,'2050-2.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2060-2.0'); n1 = find(~cellfun(@isempty,n1));
Unenth = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2070-2.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2080-2.0'); n1 = find(~cellfun(@isempty,n1));
Tired = table2array(MH(:,n1)); clear n1
RecentMDD = sum([Mood Unenth Tense Tired],2);
D4 = [Mood Unenth Tense Tired];

% PHQ
n1 = strfind(H,'20514-0.0'); n1 = find(~cellfun(@isempty,n1));
Interest = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20510-0.0'); n1 = find(~cellfun(@isempty,n1));
Down = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20517-0.0'); n1 = find(~cellfun(@isempty,n1));
Sleep = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20519-0.0'); n1 = find(~cellfun(@isempty,n1));
Tired2 = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20511-0.0'); n1 = find(~cellfun(@isempty,n1));
Appetite = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20507-0.0'); n1 = find(~cellfun(@isempty,n1));
Bad = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20508-0.0'); n1 = find(~cellfun(@isempty,n1));
Concentrate = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20518-0.0'); n1 = find(~cellfun(@isempty,n1));
Restless = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20513-0.0'); n1 = find(~cellfun(@isempty,n1));
Dead = table2array(MH(:,n1)); clear n1
PHQ = sum([Interest Down Sleep Tired2 Appetite Bad Concentrate Restless Dead],2)-9;
Dep_all = [Mood Unenth Tense Tired Interest Down Sleep Tired2 Appetite Bad Concentrate Restless Dead];

% Neuroticism
n1 = strfind(H,'1920-2.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1930-2.0'); n1 = find(~cellfun(@isempty,n1));
Miserable = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1940-2.0'); n1 = find(~cellfun(@isempty,n1));
Irritable = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1950-2.0'); n1 = find(~cellfun(@isempty,n1));
Sensitive = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1960-2.0'); n1 = find(~cellfun(@isempty,n1));
Fedup = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1970-2.0'); n1 = find(~cellfun(@isempty,n1));
Nervous = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1980-2.0'); n1 = find(~cellfun(@isempty,n1));
Worry = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'1990-2.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2000-2.0'); n1 = find(~cellfun(@isempty,n1));
Embarrassment = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2010-2.0'); n1 = find(~cellfun(@isempty,n1));
Nerves = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2020-2.0'); n1 = find(~cellfun(@isempty,n1));
Loneliness = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2030-2.0'); n1 = find(~cellfun(@isempty,n1));
Guilt = table2array(MH(:,n1)); clear n1
N = sum([Mood Miserable Irritable Sensitive Fedup Nervous Worry Tense Embarrassment Nerves Loneliness Guilt],2);

% GAD
n1 = strfind(H,'20506-0.0'); n1 = find(~cellfun(@isempty,n1));
Nervous = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20509-0.0'); n1 = find(~cellfun(@isempty,n1));
Control = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20520-0.0'); n1 = find(~cellfun(@isempty,n1));
Worry = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20515-0.0'); n1 = find(~cellfun(@isempty,n1));
Relax = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20516-0.0'); n1 = find(~cellfun(@isempty,n1));
Restless = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20505-0.0'); n1 = find(~cellfun(@isempty,n1));
Annoyed = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'20512-0.0'); n1 = find(~cellfun(@isempty,n1));
Afraid = table2array(MH(:,n1)); clear n1
GAD = sum([Nervous Control Worry Relax Restless Annoyed Afraid],2)-7;

% Combine all summary measures
Dsum = [RecurrentMDD RecentMDD PHQ N GAD]';
Dsum = nets_demean(Dsum,2);
Dsum = nets_demean((Dsum'-conf2*(Pconf2*Dsum'))',2);
Dsum = Dsum';

MH = load('CCA_inputs_behavior_CCA.mat');
CCA1 = load('CCA_results_Dsum.mat');
UV = mean([CCA1.U(:,1) CCA1.V(:,1)],2);
R = corr(MH.Dsum,UV);  R = R*-1;
figure; barh(R);
set(gca,'ytick',1:5,'yticklabel',{'Depression status','RDS_4','PHQ_9','N_1_2','GAD_7'});
title('Mental health CCA correlations');
xlabel('Cross-subject correlation with CCA score (UV)');


%% Multiply by CCA weights
CCA = load('CCA_results_DsumNEW.mat');
Vrepl = [REST_OUT IDP_OUT_NEW] * CCA.A;
Urepl = Dsum * CCA.B;
corr(Vrepl,Urepl)

%% Replilcate IDP correrllations
load('PostCCA.mat');
UV2 = mean([Urepl(:,1) Vrepl(:,1)],2);
R2 = corr([RecurrentMDD RecentMDD PHQ N GAD],UV2);
r = corr(varsd,UV2(:,1),'rows','pairwise');
r = r(Iidp); R_IDP(:,3) = r; clear r
r = corr(IDP2',UV2(:,1)); R_REST(:,3) = r; clear r

figure
spider_plot([R2 R]',...
    'AxesLabels', {'Depression status','RDS-4','PHQ-9','N-12','GAD-7'},...
    'AxesInterval', 4,...
    'AxesPrecision', 2,...
    'AxesLimits', [0,0,0,0,0;1,1,1,1,1],...
    'Color', [0, 0, 1; 1, 0, 0],...
    'AxesFontSize', 8,...
    'LabelFontSize', 10);
legend({'Correlation with UV in replication sample','Correlation with UV in exploratory CCA sample'},'Location','SouthEast')
title('B. Mental health CCA correlations');
print(gcf,'Fig_CCA_mh','-dpng','-r300');

% Find edge numbers
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
R_REST = R_REST(sig_rest,:); [~,irest] = sort(R_REST(:,1)); R_REST = R_REST(irest,:);
R_IDP = R_IDP(sig_nonrest,:); [~,inonrest] = sort(R_IDP(:,1)); R_IDP = R_IDP(inonrest,:);

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

