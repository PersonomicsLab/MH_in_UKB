clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSL')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'Subjects_CCA.csv';

% Load data
S = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
%Maps = {'SCA_amygdala','SCA_pcc','SCA_dlpfc','SCA_insula','dr_stage2'};
Maps = {};
%ICsignal = [1  2  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];
ICsignal = [1 6 7 3]; % DMN, right central executive, lefft central executive, salience (based on visual comparison against Fig 1 in Mulders et al 2015) 
Checks = zeros(5,3); % REST, SWI, TASK, DWI, STRUCTURAL
Checks_maps = zeros(5,3); % SCA_amygdala, SCA_pcc, SCA_dlpfc, SCA_insula, dr_stage2 (21x)
variance_threshold = 0.5; % keep number of PCs that explain at least 50% of variance
NkeepMax = floor(length(S)/10/2);

% Prepare confounds
load(sprintf('%s/Confounds_Subjects_%s.mat','/scratch/janine/MentalHealthInUKB/Data/',INPUT(10:end-4)));
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
    fprintf('Processing %s maps \n',Maps{m})
    maps_subs = zeros(length(GM),length(S));
    for s = 1:length(S)
        map = read_avw(sprintf('/scratch/janine/MDD_spatial/DRmaps/%s_sub-%d.nii.gz',Maps{m},S(s)));
        if size(map,4)>1
            map = map(:,:,:,ICsignal);
            for i = 1:size(map,4)
                map1 = squeeze(map(:,:,:,i));
                %maps_subs((i-1)*length(GM)+1:i*length(GM),s) = map1(GM);
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
        if Nkeep > NkeepMax; Nkeep = NkeepMax; end
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
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
REST_OUT = V(:,1:Nkeep);
Checks(1,1) = size(IDP,2); Checks(1,2) = Nkeep; Checks(1,3) = Vexp(Nkeep);
clear Nkeep V D Vexp s IDP C inds

%% Load data
DATA = readtable('/scratch/janine/MentalHealthInUKB/Data/IDP_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S); DATA = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/scratch/janine/MentalHealthInUKB/Data/ExtractVariables/vars.mat','IDP_nonrest');

%% "impute" missing  data - actually this avoids any imputation
fprintf('Processing IDPs\n');
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
varsdCOV = zeros(size(varsd,1));
for i = 1:size(varsd,1) % estimate "pairwise" covariance, ignoring missing data
    for j = 1:size(varsd,1)
        grot = varsd([i j],:); 
        grot = cov(grot(:,sum(isnan(grot))==0)'); 
        varsdCOV(i,j) = grot(1,2);
    end
end
varsdCOV2 = nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[V,D] = eig(varsdCOV2); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Vexp = diag(D)/sum(diag(D)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
IDP_OUT = V(:,1:Nkeep);
Checks(2,1) = size(varsd,2); Checks(2,2) = Nkeep; Checks(2,3) = Vexp(Nkeep);

varsd1 = knnimpute(varsd);
varsd1 = varsd1';
C = varsd1' * varsd1;
[V1,D1] = eig(C); [~,inds] = sort(diag(D1), 'descend'); D1 = D1(inds,inds); V1 = V1(:,inds);
Vexp = diag(D1)/sum(diag(D1)); Vexp = cumsum(Vexp); Nkeep = find(Vexp>=variance_threshold,1,'first');
if Nkeep > NkeepMax; Nkeep = NkeepMax; end
IDP_OUT_NEW = V1(:,1:Nkeep);
Checks(3,1) = size(varsd,2); Checks(3,2) = Nkeep; Checks(3,3) = Vexp(Nkeep);
clear Nkeep V D Vexp varsd varsdCOV varsdCOV2 i grot grotconf j

%% Save results
sprintf('Saving results\n');
save(sprintf('CCA_inputs_brain_%s.mat',INPUT(10:end-4)),'MAPS_OUT','REST_OUT','IDP_OUT','IDP_OUT_NEW','Checks','Checks_maps');


