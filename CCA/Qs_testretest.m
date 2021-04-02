clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'Subjects_testretest.csv';

% Load data
S = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
MH1 = readtable('/scratch/janine/MentalHealthInUKB/Data/MH_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(MH1(:,1)),S); MH1 = MH1(s,:);
H1 = get_UKB_headers(MH1);
MH1 = standardizeMissing(MH1,-3); MH1 = standardizeMissing(MH1,-1); MH1 = standardizeMissing(MH1,-818);
MH2 = readtable('/scratch/janine/MentalHealthInUKB/Data/MH_scan2.tsv','FileType','text');
[~,s,~] = intersect(table2array(MH2(:,1)),S); MH2 = MH2(s,:);
H2 = get_UKB_headers(MH2);
MH2 = standardizeMissing(MH2,-3); MH2 = standardizeMissing(MH2,-1); MH2 = standardizeMissing(MH2,-818);

% Probable depression status
n1 = strfind(H1,'4598-2.0'); n1 = find(~cellfun(@isempty,n1));
EverDep = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'4631-2.0'); n1 = find(~cellfun(@isempty,n1));
EverUnenth = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'4609-2.0'); n1 = find(~cellfun(@isempty,n1));
DurDep = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'5375-2.0'); n1 = find(~cellfun(@isempty,n1));
DurUnenth = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2090-2.0'); n1 = find(~cellfun(@isempty,n1));
SeenGP = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2100-2.0'); n1 = find(~cellfun(@isempty,n1));
SeenPsych = table2array(MH1(:,n1)); clear n1
Probable_depression_status_1 = zeros(size(S));
for n = 1:length(S)
    if EverDep(n) == 1 || EverUnenth(n) == 1
        if DurDep(n)>1 || DurUnenth(n)>1
            if SeenGP(n) == 1 || SeenPsych(n) == 1
                Probable_depression_status_1(n) = 1;
            end
        end
    end
end

% RDS-4
n1 = strfind(H1,'2050-2.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2060-2.0'); n1 = find(~cellfun(@isempty,n1));
Unenth = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2070-2.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2080-2.0'); n1 = find(~cellfun(@isempty,n1));
Tired = table2array(MH1(:,n1)); clear n1
RDS4_1 = sum([Mood Unenth Tense Tired],2);

% PHQ
n1 = strfind(H1,'20514-0.0'); n1 = find(~cellfun(@isempty,n1));
Interest = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20510-0.0'); n1 = find(~cellfun(@isempty,n1));
Down = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20517-0.0'); n1 = find(~cellfun(@isempty,n1));
Sleep = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20519-0.0'); n1 = find(~cellfun(@isempty,n1));
Tired2 = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20511-0.0'); n1 = find(~cellfun(@isempty,n1));
Appetite = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20507-0.0'); n1 = find(~cellfun(@isempty,n1));
Bad = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20508-0.0'); n1 = find(~cellfun(@isempty,n1));
Concentrate = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20518-0.0'); n1 = find(~cellfun(@isempty,n1));
Restless = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20513-0.0'); n1 = find(~cellfun(@isempty,n1));
Dead = table2array(MH1(:,n1)); clear n1
PHQ9_1 = sum([Interest Down Sleep Tired2 Appetite Bad Concentrate Restless Dead],2)-9;

% Neuroticism
n1 = strfind(H1,'1920-2.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1930-2.0'); n1 = find(~cellfun(@isempty,n1));
Miserable = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1940-2.0'); n1 = find(~cellfun(@isempty,n1));
Irritable = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1950-2.0'); n1 = find(~cellfun(@isempty,n1));
Sensitive = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1960-2.0'); n1 = find(~cellfun(@isempty,n1));
Fedup = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1970-2.0'); n1 = find(~cellfun(@isempty,n1));
Nervous = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1980-2.0'); n1 = find(~cellfun(@isempty,n1));
Worry = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'1990-2.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2000-2.0'); n1 = find(~cellfun(@isempty,n1));
Embarrassment = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2010-2.0'); n1 = find(~cellfun(@isempty,n1));
Nerves = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2020-2.0'); n1 = find(~cellfun(@isempty,n1));
Loneliness = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'2030-2.0'); n1 = find(~cellfun(@isempty,n1));
Guilt = table2array(MH1(:,n1)); clear n1
N12_1 = sum([Mood Miserable Irritable Sensitive Fedup Nervous Worry Tense Embarrassment Nerves Loneliness Guilt],2);

% GAD
n1 = strfind(H1,'20506-0.0'); n1 = find(~cellfun(@isempty,n1));
Nervous = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20509-0.0'); n1 = find(~cellfun(@isempty,n1));
Control = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20520-0.0'); n1 = find(~cellfun(@isempty,n1));
Worry = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20515-0.0'); n1 = find(~cellfun(@isempty,n1));
Relax = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20516-0.0'); n1 = find(~cellfun(@isempty,n1));
Restless = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20505-0.0'); n1 = find(~cellfun(@isempty,n1));
Annoyed = table2array(MH1(:,n1)); clear n1
n1 = strfind(H1,'20512-0.0'); n1 = find(~cellfun(@isempty,n1));
Afraid = table2array(MH1(:,n1)); clear n1
GAD7_1 = sum([Nervous Control Worry Relax Restless Annoyed Afraid],2)-7;

% RDS-4 scan 2
n1 = strfind(H2,'2050-3.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2060-3.0'); n1 = find(~cellfun(@isempty,n1));
Unenth = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2070-3.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2080-3.0'); n1 = find(~cellfun(@isempty,n1));
Tired = table2array(MH2(:,n1)); clear n1
RDS4_2 = sum([Mood Unenth Tense Tired],2);

% Neuroticism
n1 = strfind(H2,'1920-3.0'); n1 = find(~cellfun(@isempty,n1));
Mood = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1930-3.0'); n1 = find(~cellfun(@isempty,n1));
Miserable = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1940-3.0'); n1 = find(~cellfun(@isempty,n1));
Irritable = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1950-3.0'); n1 = find(~cellfun(@isempty,n1));
Sensitive = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1960-3.0'); n1 = find(~cellfun(@isempty,n1));
Fedup = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1970-3.0'); n1 = find(~cellfun(@isempty,n1));
Nervous = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1980-3.0'); n1 = find(~cellfun(@isempty,n1));
Worry = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'1990-3.0'); n1 = find(~cellfun(@isempty,n1));
Tense = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2000-3.0'); n1 = find(~cellfun(@isempty,n1));
Embarrassment = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2010-3.0'); n1 = find(~cellfun(@isempty,n1));
Nerves = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2020-3.0'); n1 = find(~cellfun(@isempty,n1));
Loneliness = table2array(MH2(:,n1)); clear n1
n1 = strfind(H2,'2030-3.0'); n1 = find(~cellfun(@isempty,n1));
Guilt = table2array(MH2(:,n1)); clear n1
N12_2 = sum([Mood Miserable Irritable Sensitive Fedup Nervous Worry Tense Embarrassment Nerves Loneliness Guilt],2);

RDS4_2(isnan(N12_2)) = nan;
N12_2(isnan(RDS4_2)) = nan;

length(find(isnan(RDS4_2)==0))
length(find(isnan(N12_2)==0))

%%% Save results
sprintf('Saving results\n');
save('Qs_testretest.mat','RDS4_1','RDS4_2','PHQ9_1',...
    'N12_1','N12_2','GAD7_1','Probable_depression_status_1');











