clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSL')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')
INPUT = 'CCA_results_Dsum';

% Load CCA results
CCA = load(sprintf('%s.mat',INPUT));
UV = mean([CCA.U CCA.V],2);

% Load data
S = load('/scratch/janine/MentalHealthInUKB/SubjectSplits/Subjects_CCA.csv');
Maps = {'SCA_amygdala','SCA_pcc','SCA_dlpfc','SCA_insula','dr_stage2'};
ICsignal = [1  2  3  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];

% Prepare confounds
load('/scratch/janine/MentalHealthInUKB/Data/Confounds_Subjects_CCA_28-Dec-2020.mat');
conf = nets_demean(conf);
Pconf = pinv(conf);

% Gray matter mask
%system('module load fsl')
%system('fast -o Mask/MNIseg $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz');
%system('fslmaths Mask/MNIseg_pve_1.nii.gz -thr 0 -bin Mask/GrayMatterMask.nii.gz');
GM = read_avw('Mask/GrayMatterMask.nii.gz'); GM = find(GM==1);

%% Maps
p_maps = [];
for m = 1:length(Maps)
    fprintf('Processing %s maps \n',Maps{m})
    maps_subs = zeros(length(GM),length(S));
    for s = 1:length(S)
        map = read_avw(sprintf('/scratch/janine/MDD_spatial/DRmaps/%s_sub-%d.nii.gz',Maps{m},S(s)));
        if size(map,4)>1
            map = map(:,:,:,ICsignal);
            for i = 1:size(map,4)
                map1 = squeeze(map(:,:,:,i));
                maps_subs((i-1)*length(GM)+1:i*length(GM),s) = map1(GM);
                clear map1
            end
        else
            maps_subs(:,s) = map(GM);
        end
    end
    map1 = nets_demean(maps_subs,2);
    map1 = nets_demean((map1'-conf*(Pconf*map1'))',2);
    [R,p] = corr(map1',UV); p_maps = [p_maps; p(:)];
    OUTallR = zeros(91,109,91,size(map1,1)/length(GM));
    OUTallp = zeros(91,109,91,size(map1,1)/length(GM));
    for i = 1:size(map1,1)/length(GM)
        OUT = zeros(91,109,91); 
        OUT(GM) = R((i-1)*length(GM)+1:i*length(GM)); OUTallR(:,:,:,i) = OUT;
        OUT(GM) = p((i-1)*length(GM)+1:i*length(GM)); OUTallp(:,:,:,i) = OUT;
    end
    save_avw(OUTallR,sprintf('%s_posthoc_%s_R.nii.gz',INPUT,Maps{m}),'f',[2 2 2 2]);
    save_avw(OUTallp,sprintf('%s_posthoc_%s_p.nii.gz',INPUT,Maps{m}),'f',[2 2 2 2]);
end
clear m map maps_subs s i

%% Resting state IDPs
fprintf('Processing resting state IDPs\n');
load('/scratch/janine/MentalHealthInUKB/Data/IDPs.mat');
[~,s,~] = intersect(subs,S); IDP = IDP(s,:)';
IDP = nets_demean(IDP,2);
IDP = nets_demean((IDP'-conf*(Pconf*IDP'))',2);
[R_REST,p_REST] = corr(IDP',UV);


%% Load data
DATA = readtable('/scratch/janine/MentalHealthInUKB/Data/IDP_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S); DATA = DATA(s,:);
H = get_UKB_headers(DATA); 
load('/scratch/janine/MentalHealthInUKB/Data/ExtractVariables/vars.mat','IDP');

fprintf('Processing IDPs\n');
varsd = table2array(DATA);
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
[R_IDP,p_IDP] = corr(varsd,UV,'rows','pairwise');

%% multiple comparison correction
[pthr,pcor,padj] = fdr([p_maps; p_REST; p_IDP]);

%% Save results
sprintf('Saving results\n');
save(sprintf('%s_posthoc.mat',INPUT),'R_IDP','p_IDP','R_REST','p_REST','pthr');


