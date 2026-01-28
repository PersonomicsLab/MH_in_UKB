clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'Subjects_effectsize.csv';

% Load data
S = load(sprintf('%s/%s','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/SubjectSplits/',INPUT));
MH = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/MH_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(MH(:,1)),S); MH = MH(s,:);
H = get_UKB_headers(MH);
MH = standardizeMissing(MH,-3); MH = standardizeMissing(MH,-1); MH = standardizeMissing(MH,-818);

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s.mat','/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/',INPUT(10:end-4)));
conf = nets_demean(conf);
Pconf = pinv(conf);

% Full data 
idx = ismissing(MH(:,{'x4609_2_0','x5375_2_0'}));
MH{:,{'x4609_2_0','x5375_2_0'}}(idx) = 0;
Dall = table2array(MH)';
n1 = strfind(H,'eid'); n1 = find(~cellfun(@isempty,n1));
ID = n1; clear n1
Dall(ID,:) = [];
Dall = nets_demean(Dall,2);
Dall = nets_demean((Dall'-conf*(Pconf*Dall'))',2);
Dall = Dall';

% Probable depression status
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
MDDstatus = zeros(size(S));
for n = 1:length(S)
    if EverDep(n) == 1 || EverUnenth(n) == 1
        if DurDep(n)>1 || DurUnenth(n)>1
            if SeenGP(n) == 1 || SeenPsych(n) == 1
                MDDstatus(n) = 1;
            end
        end
    end
end

% Recent depressive symptoms (RDS)
n1 = strfind(H,'2050-2.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2060-2.0'); n1 = find(~cellfun(@isempty,n1));
Unenth = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2070-2.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH(:,n1)); clear n1
n1 = strfind(H,'2080-2.0'); n1 = find(~cellfun(@isempty,n1));
Tired = table2array(MH(:,n1)); clear n1
RDS = sum([Mood Unenth Tense Tired],2);
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
Dsum = [MDDstatus RDS PHQ N GAD]';
Dsum = nets_demean(Dsum,2);
Dsum = nets_demean((Dsum'-conf*(Pconf*Dsum'))',2);
Dsum = Dsum';

%%% Save results
sprintf('Saving results\n');
save(sprintf('CCArepliation_behavior_%s.mat',INPUT(10:end-4)),'Dall','Dsum','Dep_all','D4');











