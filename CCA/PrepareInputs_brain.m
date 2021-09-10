clear all; close all; clc

addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/FSLNets')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/FSL')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/')
INPUT = 'Subjects_CCAnew.csv';

% Load data
S = load(sprintf('%s/%s','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/SubjectSplits/',INPUT));
Checks = zeros(3,3); % REST, TASK, STRUCTURAL
variance_threshold = 50; % keep number of PCs that explain at least 50% of variance
NkeepMax = floor(length(S)/10/2);

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s.mat','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/',INPUT(10:end-4)));
conf = nets_demean(conf);
Pconf = pinv(conf);

%% Resting state IDPs
fprintf('Processing resting state IDPs\n');
IDP = load('/Users/janinebijsterbosch/Box/00_WashU/Data/UKB/IDP/IDPs_rest_40k.txt');
subs = IDP(:,1); IDP = IDP(:,2:end);
[~,s,~] = intersect(subs,S); IDP = IDP(s,:);
IDP = nets_demean(IDP);
IDP = nets_demean(IDP-conf*(Pconf*IDP));
[REST_COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(IDP);
%SCORE_PROJECT = (IDP*REST_COEFF);
Vexp = cumsum(EXPLAINED); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
REST_OUT = SCORE(:,1:Nkeep);
Checks(1,1) = size(IDP,2); Checks(1,2) = Nkeep; Checks(1,3) = Vexp(Nkeep);
Nkeep_rest = Nkeep;
clear COEFF LATENT TSQUARED EXPLAINED MU Nkeep

%% Load data
DATA = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/IDP_scan1_40k.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S); DATA = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/ExtractVariables/vars.mat','IDP_nonrest');

%% Deconfound non-rest data and impute missing  data 
fprintf('Processing non resting-state IDPs\n');
varsd = table2array(DATA);
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
varsd(:,ID) = [];
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd1 = knnimpute(varsd);

%% Dimensionality reduction task IDPs
fprintf('Dimensionality reduction task\n');
task1 = strncmp('Median BOLD',IDP_nonrest(:,2),11); task1 = find(task1==1);
task2 = strncmp('Median z-stat',IDP_nonrest(:,2),13); task2 = find(task2==1);
task3 = strncmp('90th percentile',IDP_nonrest(:,2),15); task3 = find(task3==1);
task = sort([task1; task2; task3]); clear task1 task2 task3
varsd1_task = varsd1(:,task);
[TASK_COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(varsd1_task);
Vexp = cumsum(EXPLAINED); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
TASK_OUT = SCORE(:,1:Nkeep);
Checks(2,1) = size(varsd1_task,2); Checks(2,2) = Nkeep; Checks(2,3) = Vexp(Nkeep);
Nkeep_task = Nkeep;
clear COEFF LATENT TSQUARED EXPLAINED MU Nkeep

%% Dimensionality reduction on structural IDPs
fprintf('Dimensionality reduction structural\n');
nottask = setdiff(1:size(varsd1,2),task);
varsd1_nottask = varsd1(:,nottask);
[STRUCT_COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(varsd1_nottask);
Vexp = cumsum(EXPLAINED); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
STRUCT_OUT = SCORE(:,1:Nkeep);
Checks(3,1) = size(varsd1_nottask,2); Checks(3,2) = Nkeep; Checks(3,3) = Vexp(Nkeep);
Nkeep_struct = Nkeep;
clear COEFF LATENT TSQUARED EXPLAINED MU Nkeep

%% Save results
sprintf('Saving results\n');
save(sprintf('CCA_inputs_brain_%s.mat',INPUT(10:end-4)),'REST_OUT','TASK_OUT','STRUCT_OUT','Checks');
save(sprintf('Eigs_%s.mat',INPUT(10:end-4)),'REST_COEFF','Nkeep_rest','TASK_COEFF','Nkeep_task','STRUCT_COEFF','Nkeep_struct');

