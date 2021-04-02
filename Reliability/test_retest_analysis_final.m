% This script runs test-retest analyses of the brain target IDPs from the
% UK Biobank Data, and was developed during the WAPIAW 

clc
clear all
close all

% Modality codes:
% area = 1
% volume = 2
% thick = 3
% fa = 4
% md = 5
% ta = 6
% amp = 7
% fnt = 8
% pnt = 9

% modality names in order of codes
modalities = {'Area', 'Volume', 'CT', 'FA', 'MD', 'TA', 'AMP' , 'FNT', 'PNT'};

%%%% Load data files 
% navigate to data directory
data_dir = '/scratch/janine/MentalHealthInUKB/Reliability';
cd(data_dir);

% load target sheets
brain_targs_nr = readtable('WAPIAW_brain_targets.xlsx', 'Sheet', 'non_rest');
brain_targs_rs = readtable('WAPIAW_brain_targets.xlsx', 'Sheet', 'rest');
brain_targs_lt = readtable('WAPIAW_brain_targets.xlsx', 'Sheet', 'literature');

% load IDP data tables for non resting-state
idp_scan1 = readtable('IDP_scan1.tsv', 'FileType', 'delimitedtext', 'delimiter', 'tab');
idp_scan2 = readtable('IDP_scan2.tsv', 'FileType', 'delimitedtext', 'delimiter', 'tab');

% load IDP tables for resting state
idp_rest1 = load('IDPs.mat');
idp_rest2 = load('IDPs_scan2.mat');

% load confound tables
load('Confounds_Subjects_testretest_scan1.mat');
conf1 = conf; clear conf;
load('Confounds_Subjects_testretest_scan2.mat');
conf2 = conf; clear conf;

% load test-retest subject list
subj_trt = readtable('Subjects_testretest.csv');

% select out relevant subjects from rest IDP matrices
[~, ia] = intersect(idp_rest1.subs, subj_trt.Var1);
idp_rest1.IDP = idp_rest1.IDP(ia,:);
idp_rest1.subs = idp_rest1.subs(ia);

%%%% Get variable names from target lists

% non-rest target list
prefixes = repmat('x', [height(brain_targs_nr),1]); 
suffixes_1 = repmat('_2_0', [height(brain_targs_nr),1]);
suffixes_2 = repmat('_3_0', [height(brain_targs_nr),1]);
bt_nr_id_s1 = cellstr([prefixes, num2str(brain_targs_nr.id), suffixes_1]); % non-rest target IDs scan 1
bt_nr_id_s2 = cellstr([prefixes, num2str(brain_targs_nr.id), suffixes_2]); % non-rest target IDs scan 2
bt_nr_id_name = brain_targs_nr.name;

% get modalities for each
is_area_nr = contains(bt_nr_id_name, 'Area', 'IgnoreCase', true);
is_volume_nr = contains(bt_nr_id_name, 'Volume', 'IgnoreCase', true).*2;
is_thick_nr = contains(bt_nr_id_name, 'thickness', 'IgnoreCase', true).*3;
is_fa_nr = contains(bt_nr_id_name, 'FA', 'IgnoreCase', false).*4;
is_md_nr = contains(bt_nr_id_name, 'MD', 'IgnoreCase', false).*5;
is_ta_nr = contains(bt_nr_id_name, 'activation', 'IgnoreCase', true).*6;
total_count = numel(find(is_area_nr))+numel(find(is_volume_nr))+numel(find(is_thick_nr))+numel(find(is_fa_nr))+numel(find(is_md_nr))+numel(find(is_ta_nr));
disp(['Total count for non-rest types: ' num2str(total_count) ' (expected = ' num2str(size(bt_nr_id_name,1)) ')' ]);
nr_modals = is_area_nr+is_volume_nr+is_thick_nr+is_fa_nr+is_md_nr+is_ta_nr;

% rest target list
bt_rs_id_s1 = cellstr(num2str(brain_targs_rs.IndexIn3466));
bt_rs_id_s2 = bt_rs_id_s1;
bt_rs_id_name = cell(height(brain_targs_rs),1);
for i = 1:height(brain_targs_rs)
    bt_rs_id_name{i} = [brain_targs_rs.IDPType{i}(1:end-1) '_' num2str(brain_targs_rs.NodeIfAmp(i)) '_'...
        num2str(brain_targs_rs.Edge_IfEdge(i)) '_'...
        num2str(brain_targs_rs.Node1IfEdge(i)) '_'...
        num2str(brain_targs_rs.Node2IfEdge(i))];
end
is_amp_rs = contains(brain_targs_rs.IDPType(:), 'AMP').*7;
is_fnt_rs = contains(brain_targs_rs.IDPType(:), 'FNET').*8;
is_pnt_rs = contains(brain_targs_rs.IDPType(:), 'PNET').*9;
idp_set_rs = brain_targs_rs.IDPset;
total_count = numel(find(is_amp_rs))+numel(find(is_fnt_rs))+numel(find(is_pnt_rs));
disp(['Total count for rest types: ' num2str(total_count) ' (expected = ' num2str(size(bt_rs_id_name,1)) ')' ]);
rs_modals = is_amp_rs+is_fnt_rs+is_pnt_rs;

% literature target list
prefixes = repmat('x', [height(brain_targs_lt),1]); 
suffixes_1 = repmat('_2_0', [height(brain_targs_lt),1]);
suffixes_2 = repmat('_3_0', [height(brain_targs_lt),1]);
bt_lt_id_lh_s1 = cellstr([prefixes, num2str(brain_targs_lt.left), suffixes_1]); % non-rest target IDs scan 1
bt_lt_id_lh_s2 = cellstr([prefixes, num2str(brain_targs_lt.left), suffixes_2]); % non-rest target IDs scan 2
bt_lt_id_rh_s1 = cellstr([prefixes, num2str(brain_targs_lt.right), suffixes_1]); % non-rest target IDs scan 1
bt_lt_id_rh_s2 = cellstr([prefixes, num2str(brain_targs_lt.right), suffixes_2]); % non-rest target IDs scan 2
bt_lt_id_lh_name = cellstr([repmat('LH ', [height(brain_targs_lt),1]) char(brain_targs_lt.name)]);
bt_lt_id_rh_name = cellstr([repmat('RH ', [height(brain_targs_lt),1]) char(brain_targs_lt.name)]);

% get unique literature IDPs
[bt_lt_id_s1, ia1] = unique([bt_lt_id_lh_s1; bt_lt_id_rh_s1], 'stable');
[bt_lt_id_s2, ia2] = unique([bt_lt_id_lh_s2; bt_lt_id_rh_s2], 'stable');
bt_lt_id_name = [bt_lt_id_lh_name; bt_lt_id_rh_name];
bt_lt_id_name = bt_lt_id_name(ia1);
clear bt_lt_id_lh_s1 bt_lt_id_rh_s1 bt_lt_id_lh_s2 bt_lt_id_rh_s2...
    bt_lt_id_lh_name bt_lt_id_rh_name

% get types for literature IDPs
is_volume_lt = contains(bt_lt_id_name, 'Volume', 'IgnoreCase', true).*2;
is_thick_lt = contains(bt_lt_id_name, 'thickness', 'IgnoreCase', true).*3;
is_fa_lt = contains(bt_lt_id_name, 'FA', 'IgnoreCase', true).*4;
total_count = numel(find(is_volume_lt))+numel(find(is_thick_lt))+numel(find(is_fa_lt));
disp(['Total count for literature types: ' num2str(total_count) ' (expected = ' num2str(size(bt_lt_id_name,1)) ')' ]);
lt_modals = is_volume_lt + is_thick_lt + is_fa_lt;

%%%% subset data tables to extract relevant observations and variables

% get indices of included subjects and print total N for analysis
inc_subs = find(ismember(idp_scan1.eid, subj_trt.Var1)); % find included subject ID indices
disp(['Total number of subjects = ' num2str(numel(inc_subs))]); % print total N

%%%% get IDP scan 1 data for relevant subjects and variables

% remove redundant variables from one of the sets
redundant_vars_s1 = intersect(bt_nr_id_s1, bt_lt_id_s1); % overlap between non-rest and literature IDPs
[bt_nr_id_s1, ia] = setdiff(bt_nr_id_s1, redundant_vars_s1); % remove overlapping vars from non-rest
bt_nr_id_name = bt_nr_id_name(ia); % adjust non-rest name list
nr_modals = nr_modals(ia);

% get IDP scan 2 data for relevant subjects and variables
redundant_vars_s2 = intersect(bt_nr_id_s2, bt_lt_id_s2);
bt_nr_id_s2 = bt_nr_id_s2(ia);

% extract data subset of interest
idp_scan1_inc = [idp_scan1(inc_subs, [bt_nr_id_s1; bt_lt_id_s1]), array2table(idp_rest1.IDP(:,str2num(char(bt_rs_id_s1))), 'VariableNames', bt_rs_id_s1')]; % subset to extract data
idp_scan2_inc = [idp_scan2(inc_subs, [bt_nr_id_s2; bt_lt_id_s2]), array2table(idp_rest2.IDP(:,str2num(char(bt_rs_id_s2))), 'VariableNames', bt_rs_id_s2')];
idp_meta.varID1 = [bt_nr_id_s1; bt_lt_id_s1; bt_rs_id_s1];
idp_meta.varID2 = [bt_nr_id_s2; bt_lt_id_s2; bt_rs_id_s2];
idp_meta.varname = [bt_nr_id_name; bt_lt_id_name; bt_rs_id_name];
idp_meta.IDPtype = [zeros(size(bt_nr_id_s1,1),1); ones(size(bt_lt_id_s1,1),1); idp_set_rs];
idp_meta.modality = [nr_modals; lt_modals; rs_modals];

%%%%%%%% Run statistical analyses

% add Janine code to path
addpath(genpath('/scratch/janine/MentalHealthInUKB/MatlabScripts'));

% regress confounds out of columns 
include = 1:1:height(idp_scan1_inc);
disp(['Total number of included subjects after exclusions for missing data: ' num2str(numel(include))]);

% will hold confound-regressed data
res_idp_s1 = zeros(size(idp_scan1_inc{include,:}));
res_idp_s2 = zeros(size(idp_scan2_inc{include,:}));

% iterate over columns and regress out confounds
for i = 1:width(idp_scan1_inc)
    [~,~,res_idp_s1(:,i)] = regress(nets_demean(idp_scan1_inc{include,i}), [ones(numel(include), 1), conf1(include,:)]);    
    [~,~,res_idp_s2(:,i)] = regress(nets_demean(idp_scan2_inc{include,i}), [ones(numel(include), 1), conf2(include,:)]);
end

% get raw scan1-scan2 correlation
[r_idp_s1s2_raw, p_idp_s1s2_raw] = corr(idp_scan1_inc{include,:}, idp_scan2_inc{include,:}, 'rows', 'pairwise');
r_idp_s1s2_raw = diag(r_idp_s1s2_raw);

% get adjusted correlation
[r_idp_s1s2_adj, p_idp_s1s2_adj] = corr(res_idp_s1, res_idp_s2, 'rows', 'pairwise');
r_idp_s1s2_adj = diag(r_idp_s1s2_adj);

%%%%%%% Effects of time between scans

% load confound files
conf1_all = readtable('Conf_scan1.tsv', 'FileType', 'delimitedtext', 'delimiter', 'tab');
conf2_all = readtable('Conf_scan2.tsv', 'FileType', 'delimitedtext', 'delimiter', 'tab');

% get date variables for included subjects
date_1 = conf1_all{inc_subs(include), 'x53_2_0'};
date_2 = conf2_all{inc_subs(include), 'x53_3_0'};

% also get age
age_1 = conf1_all{inc_subs(include), 'x21003_2_0'};
age_2 = conf2_all{inc_subs(include), 'x21003_3_0'};

% also get sex
sex = conf1_all{inc_subs(include), 'x31_0_0'};

% get time differences and convert to days (rather than seconds)
time_dif = etime(datevec(date_2), datevec(date_1));
time_dif = time_dif./60; time_dif = time_dif./60; time_dif = time_dif./24;

% save sample stats
sample_stats.age_1.mean = mean(age_1);
sample_stats.age_1.sd = std(age_1);
sample_stats.age_2.mean = mean(age_2);
sample_stats.age_2.sd = std(age_2);
sample_stats.sex.count = numel(find(sex==1));
sample_stats.time_dif.mean = mean(time_dif);
sample_stats.time_dif.sd = std(time_dif);

% get adjusted partial correlation
[rp_idp_s1s2_adj, pp_idp_s1s2_adj] = partialcorr(res_idp_s1, res_idp_s2, time_dif, 'rows', 'pairwise');
rp_idp_s1s2_adj = diag(rp_idp_s1s2_adj);

%%%%% Plot modality and time effects on test-retest values
scrsz = get(groot, 'ScreenSize');
f = figure('Position', [455 scrsz(4)/24 scrsz(3)/1 scrsz(4)/1.5]);

% Plot distribution of time between scans 
subplot(2,2,1);
histogram(time_dif, 'binWidth', 5);
box off; grid on; ylabel('Number of Subjects'); xlabel('Days between Scans');
title('Inter-scan Interval Distribution');

% Get means and SEMs for each modality
modal_mean_adj = zeros(length(unique(idp_meta.modality)),2);
modal_sem_adj = modal_mean_adj;
for i = 1:length(unique(idp_meta.modality))
    modal_mean_adj(i,1) = nanmean(r_idp_s1s2_adj(idp_meta.modality == i));
    modal_sem_adj(i,1) = std(r_idp_s1s2_adj(idp_meta.modality == i))./sqrt(length(r_idp_s1s2_adj(idp_meta.modality == i)));
    modal_mean_adj(i,2) = nanmean(rp_idp_s1s2_adj(idp_meta.modality == i));
    modal_sem_adj(i,2) = std(rp_idp_s1s2_adj(idp_meta.modality == i))./sqrt(length(rp_idp_s1s2_adj(idp_meta.modality == i)));
end
% prepare data for boxchart and plot
subplot(2,2,2);
r_cat = [r_idp_s1s2_adj; rp_idp_s1s2_adj]; 
r_groups = [zeros(size(r_idp_s1s2_adj,1),1); ones(size(rp_idp_s1s2_adj,1),1)];
boxchart([idp_meta.modality;idp_meta.modality], r_cat, 'GroupByColor', r_groups)
set(gca, 'XtickLabel', modalities);
ylabel('Correlation (r)');
box off; grid on; hold on;
xlabel('Modality')
legend({'Confound Regression Only', 'Time-between Scans Also Partialled Out'})
title('Test-retest Correlations for Brain IDPs');
set(gca, 'Xtick', 1:1:9);

%%%%% get correlations between inter-scan IDP differences and time between scans
[dif_cor_time, dif_p_time] = corr(res_idp_s2-res_idp_s1, time_dif, 'rows', 'pairwise');
[~, ia] = sort(abs(dif_cor_time), 'descend');
time_cor_rank_table.r_value = dif_cor_time(ia);
time_cor_rank_table.modality = modalities(idp_meta.modality(ia))';
time_cor_rank_table.IDP_name = idp_meta.varname(ia);
time_cor_rank_table = struct2table(time_cor_rank_table);
head(time_cor_rank_table, 10)

%%%%% get mental health variables
mh_1 = load(fullfile(data_dir, 'CCA_inputs_behavior_testretest_scan1.mat'));
mh_2 = load(fullfile(data_dir, 'CCA_inputs_behavior_testretest_scan2.mat'));
mdd_dif = mh_2.RecentMDD - mh_1.RecentMDD;
mhc_p = find(mdd_dif == 0);
mhc_n = find(mdd_dif ~= 0);

% Plot group sizes
subplot(2,2,3)
bar([numel(mhc_n), numel(mhc_p)])
box off
grid on
ylabel('Number of Subjects');
set(gca, 'XtickLabel', {'No Change', 'Change'})
title('Mental Health Change Groups')

% no mental health change subject adjusted correlation
[r_nomhc_idp_s1s2_adj, p_nomhc_idp_s1s2_adj] = corr(res_idp_s1(mhc_n,:), res_idp_s2(mhc_n,:), 'rows', 'pairwise');
r_nomhc_idp_s1s2_adj = diag(r_nomhc_idp_s1s2_adj);

% mental health change subject adjusted correlations
[r_mhc_idp_s1s2_adj, p_mhc_idp_s1s2_adj] = corr(res_idp_s1(mhc_p,:), res_idp_s2(mhc_p,:), 'rows', 'pairwise');
r_mhc_idp_s1s2_adj = diag(r_mhc_idp_s1s2_adj);

% Plot test-retest correlations by sub-group
subplot(2,2,4);
r_cat = [r_nomhc_idp_s1s2_adj; r_mhc_idp_s1s2_adj]; 
r_groups = [zeros(size(r_nomhc_idp_s1s2_adj,1),1); ones(size(r_mhc_idp_s1s2_adj,1),1)];
boxchart([idp_meta.modality;idp_meta.modality], r_cat, 'GroupByColor', r_groups)
set(gca, 'XtickLabel', modalities);
ylabel('Correlation (r)');
box off; grid on; hold on;
xlabel('Modality')
legend({'No Change', 'Change'})
title('Test-retest Correlations for Brain IDPs');
set(gca, 'Xtick', 1:1:9);

%%%%% get correlations between inter-scan IDP differences and MHC change
[dif_cor_mhc, dif_p_mhc] = corr(res_idp_s2-res_idp_s1, mdd_dif, 'rows', 'pairwise');
[~, ia] = sort(abs(dif_cor_mhc), 'descend');
mhc_cor_rank_table.r_value = dif_cor_mhc(ia);
mhc_cor_rank_table.modality = modalities(idp_meta.modality(ia))';
mhc_cor_rank_table.IDP_name = idp_meta.varname(ia);
mhc_cor_rank_table = struct2table(mhc_cor_rank_table);
head(mhc_cor_rank_table, 10)

idp_difs = (res_idp_s2-res_idp_s1).^2;
[h, p, ci, stats] = ttest2(idp_difs(mhc_n,:), idp_difs(mhc_p,:));