clear all; close all; clc

addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/')
warning('off')

%% Load all possible subjects
S20 = load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/ExtractVariables/subs20k.csv');
S40 = load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/ExtractVariables/subs40k.csv');
S = setdiff(S40,S20); %clear S40

%% Check that no existing subjects are included
Scca = load('Subjects_CCA.csv'); Sretest = load('Subjects_effectsize.csv');
S = setdiff(S,Scca); S = setdiff(S,Sretest);

%% Remove subjects with missing imaging data
subs = load('/Users/janinebijsterbosch/Box/00_WashU/Data/UKB/IDP/subs40k.csv');
[S,~,~] = intersect(subs,S); %clear subs
fprintf('Subjects with imaging data: %d\n',length(S))
fprintf('with previous sample Subjects with imaging data: %d\n',length(S)+14615)

%% Remove subjects with head motion > 0.2
Conf1 = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/Conf_scan1_40k.tsv','FileType','text');
Conf1H = get_UKB_headers(Conf1); 
[~,i,~] = intersect(table2array(Conf1(:,1)),S); Conf1 = Conf1(i,:); clear i
n1 = strfind(Conf1H,'25741-2.0'); n1 = find(~cellfun(@isempty,n1));
motion_r1 = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf1H,'25742-2.0'); n1 = find(~cellfun(@isempty,n1));
motion_t1 = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf1H,'31-0.0'); n1 = find(~cellfun(@isempty,n1));
sex = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf1H,'21003-2.0'); n1 = find(~cellfun(@isempty,n1));
age = table2array(Conf1(:,n1)); clear n1
Sout = unique([find(motion_r1>0.2); find(motion_t1>0.2)]);
S(Sout) = [];
Conf1(Sout,:) = []; motion_r1(Sout) = []; motion_t1(Sout) = []; age(Sout) = []; sex(Sout) = [];
fprintf('Subjects with low motion: %d (deleted due to motion > 0.2: %d)\n',length(S),length(Sout));
fprintf('with previous sample Subjects with low motion: %d (deleted due to motion > 0.2: %d)\n',length(S)+12290,length(Sout)+2325);
clear Sout

%% Remove subjects with missing mental health data
MH = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/MH_scan1_40k.tsv','FileType','text');
MHH = get_UKB_headers(MH); 
[~,i,~] = intersect(table2array(MH(:,1)),S); 
MH = MH(i,:); clear i
MH = standardizeMissing(MH,-3); MH = standardizeMissing(MH,-1); MH = standardizeMissing(MH,-818);
Sout = ismissing(MH); Sout(:,[21 23]) = 0;
Sout = sum(Sout,2);
Sout = find(Sout>0);
S(Sout) = []; 
Conf1(Sout,:) = []; motion_r1(Sout) = []; motion_t1(Sout) = []; age(Sout) = []; sex(Sout) = []; MH(Sout,:) = [];
fprintf('Subjects with MH: %d (deleted due to missing MH: %d)\n',length(S),length(Sout));
fprintf('with previous sample Subjects with MH: %d (deleted due to missing MH: %d)\n',length(S)+8771,length(Sout)+3519);
clear Sout

%% Subset of lateonset depression
E = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/Extras_40k.tsv','FileType','text');
EH = get_UKB_headers(E); 
[~,i,~] = intersect(table2array(E(:,1)),S); 
E = E(i,:); clear i
E = standardizeMissing(E,-121); E = standardizeMissing(E,-818);
n1 = strfind(EH,'20433-0.0'); n1 = find(~cellfun(@isempty,n1));
AgeOnset = table2array(E(:,n1)); clear n1
Sout = find(AgeOnset>=60); Slateonset = S(Sout);
motionLOD = motion_r1(Sout); sexLOD = sex(Sout); ageLOD = age(Sout);
S(Sout) = [];
Conf1(Sout,:) = []; motion_r1(Sout) = []; motion_t1(Sout) = []; age(Sout) = []; sex(Sout) = [];  MH(Sout,:) = [];
fprintf('Subjects with normal onset depression: %d (subs used for late onset: %d)\n',length(S),length(Sout));
fprintf('with previous sample Subjects with normal onset depression: %d (subs used for late onset: %d)\n',length(S)+7921,length(Sout)+226);
clear Sout

%% Find matched HC group for MH group
n1 = strfind(MHH,'2090-2.0'); n1 = find(~cellfun(@isempty,n1));
SeenGP = table2array(MH(:,n1)); clear n1

subsSeenGP = find(SeenGP==1);
motion0a = motion_r1(subsSeenGP); sex0a = sex(subsSeenGP); age0a = age(subsSeenGP);
subsNoGP = find(SeenGP==0);
motion0b = motion_r1(subsNoGP); sex0b = sex(subsNoGP); age0b = age(subsNoGP);
subs0b_new = zeros(size(subsSeenGP)); 
for n = 1:length(subsSeenGP)
    s1 = find(sex0b==sex0a(n));
    agediff = abs(age0b(s1)-age0a(n)); 
    s2 = find(agediff==min(agediff)); s1=s1(s2);
    [~,s] = min(abs((motion0b(s1)-motion0a(n)))); s = s1(s);
    subs0b_new(n) = subsNoGP(s); 
    sex0b(s) = nan; age0b(s) = nan; motion0b(s) = nan;
    clear s s1 s2 agediff;
end

% Match late onset subjects
subsLOD_HC = zeros(size(Slateonset)); 
for n = 1:length(Slateonset)
    s1 = find(sex0b==sexLOD(n));
    agediff = abs(age0b(s1)-ageLOD(n)); 
    s2 = find(agediff==min(agediff)); s1=s1(s2);
    [~,s] = min(abs((motion0b(s1)-motionLOD(n)))); s = s1(s);
    subsLOD_HC(n) = subsNoGP(s); 
    sex0b(s) = nan; age0b(s) = nan; motion0b(s) = nan;
    clear s s1 s2 agediff;
end
subsSeenGP = S(subsSeenGP); 
subsSeenGP_controls = S(subs0b_new); 
Slateonset_controls = S(subsLOD_HC);
motion0b = motion_r1(subs0b_new); sex0b = sex(subs0b_new); age0b = age(subs0b_new);
clear subs0b_new subs_LOD_HC;
%[H,P,CI,STATS] = ttest2(age0a,age0b)
%[H,P,CI,STATS] = ttest2(motion0a,motion0b)
fprintf('Subjects who saw GP for anxiety/depression: %d\n',length(subsSeenGP));
fprintf('with previous sample Subjects who saw GP for anxiety/depression: %d\n',length(subsSeenGP)+2426);
fprintf('Discarded subjects who never saw GP for anxiety/depression: %d\n',length(subsNoGP)-length(Slateonset_controls)-length(subsSeenGP_controls));
fprintf('with previous sample Discarded subjects who never saw GP for anxiety/depression: %d\n',length(subsNoGP)-length(Slateonset_controls)-length(subsSeenGP_controls)+2843);

Scca_new = [subsSeenGP; subsSeenGP_controls];
Scca_original = load('Subjects_cca.csv');
Scca = sort([Scca_original; Scca_new]);

%% Write out csv files 
%dlmwrite('Subjects_CCAnew.csv',Scca,'precision',7);
warning('on')



