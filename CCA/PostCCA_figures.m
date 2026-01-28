clear all; close all; clc

addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/FSLNets')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/FSL')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/')
addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/FastICA_25')
INPUT = 'CCA_results_CCAnew'; 

% Load CCA results
CCA = load(sprintf('%s.mat',INPUT)); 
UV = mean([CCA.U(:,1) CCA.V(:,1)],2);
S = load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/SubjectSplits/Subjects_CCAnew.csv');

% Prepare confounds
load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/Confounds_Subjects_CCAnew.mat');
conf = nets_demean(conf);
Pconf = pinv(conf);

% Load IDP name and order information
load('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/ExtractVariables/vars.mat')
IDPnames = IDP_nonrest;
load('CCA_inputs_behavior_CCAnew.mat');

% Reorder mental health
for n = 1:38; MHold(n) = str2double(MH(n,1)); end
MHnew = sort(MHold);
for n = 1:38; i = find(MHnew==MHold(n)); Imh(n) = i; end

% Load non-rest IDPs
DATA = readtable('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/Data/IDP_scan1_40k.tsv','FileType','text');
[~,s,~] = intersect(table2array(DATA(:,1)),S); DATA = DATA(s,:);
H = get_UKB_headers(DATA); 

% Reorder non-rest IDPs
H = H(2:end);
for n = 1:length(IDPnames); i = find(strncmp(IDPnames(n,1),H,5)==1); Iidp(n) = i; end
clear n i

% deconfound
varsd = table2array(DATA);
varsd = palm_inormal(varsd); % Gaussianise
for i = 1:size(varsd,2) % deconfound ignoring missing data
    grot = (isnan(varsd(:,i))==0); 
    grotconf = nets_demean(conf(grot,:)); 
    varsd(grot,i) = normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
clear grot grontconf i

% Load Resting state IDPs
fprintf('Processing resting state IDPs\n');
IDP = load('/Users/janinebijsterbosch/Box/00_WashU/Data/UKB/IDP/IDPs_rest_40k.txt');
subs = IDP(:,1); IDP = IDP(:,2:end);
[~,s,~] = intersect(subs,S); IDP = IDP(s,:);
IDP = nets_demean(IDP);
IDP = nets_demean(IDP-conf*(Pconf*IDP));

for n = 1
    % Correlation againts mental health
    figure; set(gcf,'Position',[100 400 1800 600],'PaperPositionMode','auto');
    [R,p] = corr(UV(:,n),Dall); R = R(Imh); p = p(Imh);
    subplot(1,2,1); barh(R);
    set(gca,'ytick',1:38,'yticklabel',MH(:,2));
    [R,p] = corr(UV(:,n),Dsum);
    subplot(1,2,2); barh(R);
    set(gca,'ytick',1:5,'yticklabel',{'RecurrentMDD', 'RecentMDD', 'PHQ', 'N', 'GAD'})
    
    [R_IDP,p_IDP] = corr(varsd,UV(:,n),'rows','pairwise');
    R_IDP = R_IDP(Iidp); p_IDP = p_IDP(Iidp);
    [R_REST,p_REST] = corr(IDP,UV(:,n));
    %[pthr,pcor,padj] = fdr([p_REST; p_IDP]);
    pthr = 0.05/length([p_REST; p_IDP]); % Bonferroni
    
    % Colors for non-rests IDPs
    C(1:153) = 1; % Volumes
    C(154:215) = 2; % Area
    C(216:277) = 3; % Mean thickness
    C(278) = 4; % hyperintensities
    C(279:305) = 5; % Weighted-mean FA
    C(306:332) = 6; % Weighted-mean MD
    C(333:348) = 7; % Task
    C(349:362) = 8; % Median T2star
    
    figure; set(gcf,'Position',[1 300 1500 800],'PaperPositionMode','auto');
    subplot('position',[0.3 0.4 0.65 0.55])
    scatter(1:length(R_IDP(:,1)),R_IDP(:,1),10,C,'filled')
    title('CCA results for non-rest IDPs');
    %ylabel('Correlation between mean UV and IDP');
    I = find(p_IDP<pthr); sig_nonrest = I;
    set(gca,'xtick',I,'xticklabel',IDP_nonrest(I,2)); XYrotalabel;
    
    % Colors for rest IDPs
    C(1:55) = 1; % 100 amps
    C(56:1540) = 2; % 100 Fnets
    C(1541:3025) = 3; % 100 Pnets
    C(3026:3046) = 4; % 25 amps
    C(3047:3256) = 5; % 25 Fnets
    C(3257:3466) = 6; % 25 Pnets
    
    figure; set(gcf,'Position',[500 300 1500 800],'PaperPositionMode','auto');
    subplot('position',[0.3 0.4 0.65 0.55])
    scatter(1:length(R_REST(:,1)),R_REST(:,1),10,C,'filled')
    title('CCA results for rest IDPs');
    %ylabel('Correlation between mean UV and IDP');
    I = find(p_REST<pthr); 
    %I = find(abs(R_REST(:,1))>0.1);
    sig_rest = I;
    %set(gca,'xtick',I,'xticklabel',IDP_rest(I,:)); XYrotalabel;
    clear C
end
sig_nonrest_names = IDP_nonrest(sig_nonrest,2);
sig_nonrest_varID = IDP_nonrest(sig_nonrest,1);
sig_rest_names = IDP_rest(sig_rest,:);
%save('PostCCA.mat','sig_nonrest','sig_rest','Iidp','sig_nonrest_names','sig_nonrest_varID','sig_rest_names','R_IDP','R_REST');

