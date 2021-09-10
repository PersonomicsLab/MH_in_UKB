clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSL')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'Subjects_CCAlatelife.csv';

% Load data
S = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
Maps = {'SCA_amygdala','SCA_pcc','SCA_dlpfc','SCA_insula','dr_stage2'};
ICsignal = [1  2  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];
Checks = zeros(5,3); % REST, SWI, TASK, DWI, STRUCTURAL
Checks_maps = zeros(25,3); % SCA_amygdala, SCA_pcc, SCA_dlpfc, SCA_insula, dr_stage2 (21x)
variance_threshold = 0.25; % keep number of PCs that explain at least 50% of variance

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s_%s.mat','/scratch/janine/MentalHealthInUKB/Data/',INPUT(10:end-4),date));
conf = nets_demean(conf);
Pconf = pinv(conf);

% Gray matter mask
%system('module load fsl')
%system('fast -o Mask/MNIseg $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz');
%system('fslmaths Mask/MNIseg_pve_1.nii.gz -thr 0 -bin Mask/GrayMatterMask.nii.gz');
GM = read_avw('Mask/GrayMatterMask.nii.gz'); GM = find(GM==1);

%% Maps
MAPS_OUT = []; I = 1;
for m = 1:length(Maps)
    fprintf('Processing %s maps \n',Maps{m});
    if strcmp(Maps{m},'dr_stage2')
        maps_subs = zeros(length(GM),length(S),length(ICsignal));
    else
        maps_subs = zeros(length(GM),length(S));
    end
    for s = 1:length(S)
        map = read_avw(sprintf('/scratch/janine/MDD_spatial/DRmaps/%s_sub-%d.nii.gz',Maps{m},S(s)));
        if size(map,4)>1
            map = map(:,:,:,ICsignal);
            for i = 1:size(map,4)
                map1 = squeeze(map(:,:,:,i));
                maps_subs(:,s,i) = map1(GM);
                clear map1
            end
        else
            maps_subs(:,s) = map(GM);
        end
    end
    for i = 1:size(maps_subs,3)
        map1 = squeeze(maps_subs(:,:,i));
        map1 = nets_demean(map1,2);
        map1 = nets_demean((map1'-conf*(Pconf*map1'))',2);
        C = map1' * map1;
        [V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
        Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
        Checks_maps(I,1) = size(GM,1); Checks_maps(I,2) = Nkeep; Checks_maps(I,3) = Vexp(Nkeep); I = I+1;
        MAPS_OUT = [MAPS_OUT V(:,1:Nkeep)];        
        clear map1 C V D inds Vexp Nkeep
    end
end
clear m map maps_subs s i

%% Resting state IDPs
fprintf('Processing resting state IDPs\n');
load('/scratch/janine/MentalHealthInUKB/Data/IDPs.mat');
[~,s,~] = intersect(subs,S); IDP = IDP(s,:)';
IDP = nets_demean(IDP,2);
IDP = nets_demean((IDP'-conf*(Pconf*IDP'))',2);
C = IDP' * IDP;
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
REST_OUT = V(:,1:Nkeep);
Checks(1,1) = size(IDP,2); Checks(1,2) = Nkeep; Checks(1,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s IDP C inds

%% SWI IDPs
fprintf('Processing SWI IDPs\n');
DATA = readtable('/scratch/janine/MentalHealthInUKB/Data/IDP_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S); DATA = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/scratch/janine/MentalHealthInUKB/Data/ExtractVariables/vars.mat','IDP');
I = strfind(IDP(:,2),'T2star');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
nall = [];
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
SWI = table2array(DATA(:,nall))'; 
SWI = nets_demean(SWI,2);
SWI = nets_demean((SWI'-conf*(Pconf*SWI'))',2);
C = SWI' * SWI;
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
SWI_OUT = V(:,1:Nkeep);
Checks(2,1) = size(SWI,2); Checks(2,2) = Nkeep; Checks(2,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s C inds n n1 nall I 

%% TASK IDPs
sprintf('Processing Task IDPs\n');
I = strfind(IDP(:,2),'BOLD effect');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
nall = [];
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'z-statistic');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
TASK = table2array(DATA(:,nall))'; 
TASK = nets_demean(TASK,2);
TASK = nets_demean((TASK'-conf*(Pconf*TASK'))',2);
C = TASK' * TASK;
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
TASK_OUT = V(:,1:Nkeep);
Checks(3,1) = size(TASK,2); Checks(3,2) = Nkeep; Checks(3,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s C inds n n1 nall I 

%% DWI IDPs
fprintf('Processing DWI IDPs\n');
I = strfind(IDP(:,2),'ean FA');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
nall = [];
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'ean ICVF');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'ean ISOVF');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'ean L');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'ean MD');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'ean OD');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
DWI = table2array(DATA(:,nall))'; 
DWI = nets_demean(DWI,2);
DWI = nets_demean((DWI'-conf*(Pconf*DWI'))',2);
C = DWI' * DWI;
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
DWI_OUT = V(:,1:Nkeep);
Checks(4,1) = size(DWI,2); Checks(4,2) = Nkeep; Checks(4,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s C inds n n1 nall I 

%% STRUCTURAL IDPs
fprintf('Processing structural IDPs\n');
I = strfind(IDP(:,2),'olume');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
nall = [];
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'thickness');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
I = strfind(IDP(:,2),'Area');
I = find(~cellfun(@isempty,I));
I = IDP(I,1);
for n = 1:length(I)
    n1 = strfind(H,sprintf('%s-2.0',I(n))); 
    nall = [nall; find(~cellfun(@isempty,n1))]; clear n1
end
STRUCT = table2array(DATA(:,nall))'; 
STRUCT = nets_demean(STRUCT,2);
STRUCT = nets_demean((STRUCT'-conf*(Pconf*STRUCT'))',2);
C = STRUCT' * STRUCT;
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
STRUCT_OUT = V(:,1:Nkeep);
Checks(5,1) = size(STRUCT,2); Checks(5,2) = Nkeep; Checks(5,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s C inds n n1 nall I 

%% Save results
sprintf('Saving results\n');
save(sprintf('CCA_inputs_brain_%s_%s.mat',INPUT(10:end-4),date),'MAPS_OUT','REST_OUT','SWI_OUT','TASK_OUT','DWI_OUT','STRUCT_OUT','Checks','Checks_maps');


