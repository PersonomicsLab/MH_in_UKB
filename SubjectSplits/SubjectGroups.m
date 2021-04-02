clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
warning('off')

%% Load all possible subjects
Sall = load('/scratch/janine/MentalHealthInUKB/Data/ExtractVariables/subs20k.csv');

%% Remove subjects with missing imaging data
I = dir('/scratch/janine/MDD_spatial/DRmaps/dr_stage2_sub-*.nii.gz');
Sica = zeros(length(I),1);
for n = 1:length(I)
    Sica(n) = str2double(I(n).name(15:21));
end
clear I n
I = dir('/scratch/janine/MDD_spatial/DRmaps/SCA_pcc_sub-*.nii.gz');
Sseed = zeros(length(I),1);
for n = 1:length(I)
    Sseed(n) = str2double(I(n).name(13:19));
end
clear I n
S = intersect(Sica,Sseed);
S = intersect(Sall,S);
clear Sall Sica Sseed
fprintf('Subjects with imaging data: %d\n',length(S))

%% Remove subjects with head motion > 0.2
Conf1 = readtable('/scratch/janine/MentalHealthInUKB/Data/Conf_scan1.tsv','FileType','text');
Conf1H = get_UKB_headers(Conf1); 
[~,i,~] = intersect(table2array(Conf1(:,1)),S); Conf1 = Conf1(i,:); clear i
Conf2 = readtable('/scratch/janine/MentalHealthInUKB/Data/Conf_scan2.tsv','FileType','text');
Conf2H = get_UKB_headers(Conf2); 
[~,i,~] = intersect(table2array(Conf2(:,1)),S); Conf2 = Conf2(i,:); clear i
n1 = strfind(Conf1H,'25741-2.0'); n1 = find(~cellfun(@isempty,n1));
motion_r1 = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf2H,'25741-3.0'); n1 = find(~cellfun(@isempty,n1));
motion_r2 = table2array(Conf2(:,n1)); clear n1
n1 = strfind(Conf1H,'25742-2.0'); n1 = find(~cellfun(@isempty,n1));
motion_t1 = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf2H,'25742-3.0'); n1 = find(~cellfun(@isempty,n1));
motion_t2 = table2array(Conf2(:,n1)); clear n1
n1 = strfind(Conf1H,'31-0.0'); n1 = find(~cellfun(@isempty,n1));
sex = table2array(Conf1(:,n1)); clear n1
n1 = strfind(Conf1H,'21003-2.0'); n1 = find(~cellfun(@isempty,n1));
age = table2array(Conf1(:,n1)); clear n1
Sout = unique([find(motion_r1>0.2); find(motion_r2>0.2); find(motion_t1>0.2); find(motion_t2>0.2)]);
S(Sout) = [];
Conf1(Sout,:) = []; Conf2(Sout,:) = []; motion_r1(Sout) = []; motion_r2(Sout) = []; motion_t1(Sout) = []; motion_t2(Sout) = []; age(Sout) = []; sex(Sout) = [];
fprintf('Subjects with low motion: %d (deleted due to motion > 0.2: %d)\n',length(S),length(Sout));
clear Sout

%% Remove subjects with missing mental health data
MH = readtable('/scratch/janine/MentalHealthInUKB/Data/MH_scan1.tsv','FileType','text');
MHH = get_UKB_headers(MH); 
[~,i,~] = intersect(table2array(MH(:,1)),S); 
MH = MH(i,:); clear i
MH = standardizeMissing(MH,-3); MH = standardizeMissing(MH,-1); MH = standardizeMissing(MH,-818);
Sout = ismissing(MH); Sout(:,[21 23]) = 0;
Sout = sum(Sout,2);
Sout = find(Sout>0);
S(Sout) = []; 
Conf1(Sout,:) = []; Conf2(Sout,:) = []; motion_r1(Sout) = []; motion_r2(Sout) = []; motion_t1(Sout) = []; motion_t2(Sout) = []; age(Sout) = []; sex(Sout) = []; MH(Sout,:) = [];
fprintf('Subjects with MH: %d (deleted due to missing MH: %d)\n',length(S),length(Sout));
clear Sout

%% Select test-retest subjects
I = dir('/scratch/janine/MDD_spatial/DRmaps/scan2_dr_stage2_sub-*.nii.gz');
Sica = zeros(length(I),1);
for n = 1:length(I)
    Sica(n) = str2double(I(n).name(21:27));
end
clear I n
I = dir('/scratch/janine/MDD_spatial/DRmaps/scan2_SCA_pcc_sub-*.nii.gz');
Sseed = zeros(length(I),1);
for n = 1:length(I)
    Sseed(n) = str2double(I(n).name(19:25));
end
clear I n
Sreliability = intersect(Sica,Sseed);
clear Sica Sseed Sout
[Sreliability,Sout,~] = intersect(S,Sreliability);
S(Sout) = [];
Conf1(Sout,:) = []; Conf2(Sout,:) = []; motion_r1(Sout) = []; motion_r2(Sout) = []; motion_t1(Sout) = []; motion_t2(Sout) = []; age(Sout) = []; sex(Sout) = [];  MH(Sout,:) = [];
fprintf('Subjects without second scan: %d (subs used for relatibility: %d)\n',length(S),length(Sout));
clear Sout

%% Subset of lateonset depression
E = readtable('/scratch/janine/MentalHealthInUKB/Data/Extras.tsv','FileType','text');
EH = get_UKB_headers(E); 
[~,i,~] = intersect(table2array(E(:,1)),S); 
E = E(i,:); clear i
E = standardizeMissing(E,-121); E = standardizeMissing(E,-818);
n1 = strfind(EH,'20433-0.0'); n1 = find(~cellfun(@isempty,n1));
AgeOnset = table2array(E(:,n1)); clear n1
Sout = find(AgeOnset>=60); Slateonset = S(Sout);
motionLOD = motion_r1(Sout); sexLOD = sex(Sout); ageLOD = age(Sout);
S(Sout) = [];
Conf1(Sout,:) = []; Conf2(Sout,:) = []; motion_r1(Sout) = []; motion_r2(Sout) = []; motion_t1(Sout) = []; motion_t2(Sout) = []; age(Sout) = []; sex(Sout) = [];  MH(Sout,:) = [];
fprintf('Subjects with normal onset depression: %d (subs used for late onset: %d)\n',length(S),length(Sout));
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
fprintf('Discarded subjects who never saw GP for anxiety/depression: %d\n',length(subsNoGP)-length(Slateonset_controls)-length(subsSeenGP_controls));


%% Split sample into exploratory CCA and confirmatory effect sizes
I = randperm(length(subsSeenGP));
subsSeenGP = subsSeenGP(I); subsSeenGP_controls = subsSeenGP_controls(I);
N = floor(length(subsSeenGP)/2);
Scca = [subsSeenGP(1:N); subsSeenGP_controls(1:N)];
Seffect =[subsSeenGP(end-N+1:end); subsSeenGP_controls(end-N+1:end)];

%% Write out csv files 
dlmwrite('Subjects_CCA.csv',Scca,'precision',7);
dlmwrite('Subjects_CCApatients.csv',subsSeenGP(end-N+1:end),'precision',7);
ageCCA = zeros(size(Scca)); for n = 1:length(Scca); ageCCA(n) = age(S==Scca(n)); end 
Slld = Scca(ageCCA>=60);
dlmwrite('Subjects_CCAlatelife.csv',Slld,'precision',7);
dlmwrite('Subjects_CCAlatelifepatient.csv',Slld(1:size(Slld,1)/2),'precision',7);
dlmwrite('Subjects_effectsize.csv',Seffect,'precision',7);
dlmwrite('Subjects_testretest.csv',Sreliability,'precision',7);
dlmwrite('Subjects_lateonset.csv',[Slateonset Slateonset_controls],'precision',7);
warning('on')



