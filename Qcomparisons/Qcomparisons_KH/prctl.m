% % %%after preparing behavioral inputs, run percentiles on them
% % Dep_all_prctl = prctile(Dep_all,[0,100]);
% % N_prctl = prctile(N, [0,100]);
% % GAD_prctl = prctile(GAD, [0,100]);
% % Dep_sum_prctl = prctile(Dep_sum, [0,100]);
% % 
% % 
% % %%scatter plot and fit line 
% % %scatter(D4, Dep_all)
% % RecentMDD_prctl = prctile(RecentMDD,[1:100]);
% % %RecentMDD_prctl = prctile(RecentMDD,[0,100], all)
% % %RecentMDD_prctl = prctile(RecentMDD,[0,100], 2)
% % %RecentMDD_prctl = prctile(RecentMDD,5, 2)
% % %RecentMDD_prctl = prctile(RecentMDD,5)
% % %scatter(RecentMDD, PHQ)
% % %scatter(PHQ, RecentMDD)
% % %PHQ_RMMDlm = fitlm(PHQ, RecentMDD)

%% histograms
histogram(PHQ)
title('PHQ')
figure, histogram(RecentMDD)
title('RDS')
figure, histogram(N)
title('Neuroticism')
figure, histogram(GAD)
title('GAD')

%% PHQ v RDS
% %RDS unique quantiles applied to PHQ (since RDS is not as granular - it's
% %the rate limiting factor if you will
%%PHQ has 18 options, RDS has 10? So RDS cut offs goes into PHQ
RDS_prctl = prctile(RecentMDD,[1:100]);
PHQ_prctl = prctile(PHQ,[1:100]);
[PHQ_unq_ptl,PHQ_cutoffs] = prctile_unique(PHQ);
[RDS_unq_ptl,RDS_cutoffs] = prctile_unique(RecentMDD);
PHQ_cutoffsrds = prctile(PHQ, RDS_unq_ptl);

%% Graph PHQ v RDS
figure, subplot(1,2,1)
scatter(PHQ_prctl, RDS_prctl)
xlabel('PHQ')
ylabel('RDS')
title('PHQ vs RDS (percentiles)')

subplot(1,2,2)
scatter(PHQ_cutoffsrds, RDS_cutoffs)
xlabel('PHQ (by RDS cutoff)')
ylabel('RDS (by RDS cutoff)')
title('PHQ vs RDS (by RDS cutoff quantiles)')

%% GAD v Neuroticism
%%GAD has 15 options, N has 13. So N cutoffs goes into GAD
GAD_prctl = prctile(GAD,[1:100]);
N_prctl = prctile(N,[1:100]);

[N_unq_ptl,N_cutoffs] = prctile_unique(N);
[GAD_unq_ptl,GAD_cutoffs] = prctile_unique(GAD);
GAD_cutoffsN = prctile(GAD, N_unq_ptl);

%for funsies:
%N_cutoffsgad = prctile(N, GAD_unq_ptl);
%scatter(N_cutoffsgad, GAD_cutoffs)

%% Graph GAD v Neuroticism
figure, subplot(1,2,1)
scatter(GAD_prctl, N_prctl)
xlabel('GAD')
ylabel('Neuroticism')
title('GAD vs Neuroticism (percentiles)')

subplot(1,2,2)
scatter(GAD_cutoffsN, N_cutoffs)
xlabel('GAD (by N cutoff)')
ylabel('Neuroticism (by N cutoff)')
title('GAD vs Neuroticism (by N cutoff quantiles)')

%%
% %Neuroticism cut-off equivalent of GAD 9
% %okay so GAD9 is the 94th percentile is between 10 and 11 in neuroticism
% %GAD9 = N10-11
% %HOW DO I GRAPH THIS
% 
% %RDS cut-off equivalent of PHQ-6 and of PHQ-15 
% %PHQ6 = RDS8-9, PHQ15 = RDS14
% 
% %PHQ6 is row/column 8 (89th percentile), PHQ15 is row/column 14 (98th percentile)
% %PHQ6 = row 10-11 in N (N=9-10),  PHQ15 = row 13 (N=12)
% %so PHQ6 = N9-10, PHQ15=N12


% %GAD9 = 95th percentile
%N95 = prctile(N,95)
% %GAD9 = N10

% %PHQ6 = 84nd percentile?
% %PHQ15 = 98th percentile?
% RDS84 = prctile(RecentMDD, 84)
% RDS98 = prctile(RecentMDD, 98)
% %PHQ6 = RDS7
% %PHQ15 = RDS12


%% RDS v Neuroticism
%RDS has 10, N has 13. So RDS cutoffs goes into N
N_cutoffsrds = prctile(N, RDS_unq_ptl);

%% Graph PHQ v RDS
figure, subplot(1,2,1)
scatter(N_prctl, RDS_prctl)
xlabel('Neuroticism')
ylabel('RDS')
title('Neuroticism vs RDS (percentiles)')

subplot(1,2,2)
scatter(N_cutoffsrds, RDS_cutoffs)
xlabel('Neuroticism (by RDS cutoff)')
ylabel('RDS (by RDS cutoff)')
title('Neuroticism vs RDS (by RDS cutoff quantiles)')


%% test retest
%righ now I have scan 2 stuff loaded so I need to change D4 to D4_2
%D4_2 = D4;
%now I've added scan 1 D4 (called D4_1)
%okay but 80% are empty in scan2 so I have to index out all the NaNs
% S2 = sum(D4_2,2);
% bad = isnan(S2);
% 
% D41_new = D4_1(~bad,:);
% D42_new = D4_2(~bad,:);
% 
% 
% [tstrtstcorr, tstrtstp] = corr(D41_new, D42_new);
% 
% S1 = sum(D4_1,2);
% RecentMDD1_new = S1(~bad, :);
% RecentMDD2_new = S2(~bad, :);
sum_corr = corr(RecentMDD1_new,RecentMDD2_new)

%% cronbach
%D41_alpha = cronbach(D4_1)
%0.8027
PHQ_a = cronbach(Dep_all)

%%
%test retest neuroticism
%N2 = N;
%then clear scan 2 and load scan 1
%N1 = N;
%then clear all but N1 and N2
% bad = isnan(N2);
% N1_sub = N1(~bad, :);
% N2_sub = N2(~bad, :);
%N12_corr = corr(N1_sub, N2_sub)


%% final figure
figure
sgtitle('Mental Health Questionnaire Comparison')
subplot(4, 2, 1)
histogram(PHQ)
title('Patient Health Questionnaire (PHQ)')
xlabel('Score')
ylabel('# of people')
subplot(4, 2, 2)
histogram(RecentMDD)
title('Recent Depression Score (RDS)')
xlabel('Score')
ylabel('# of people')
subplot(4, 2, 3)
histogram(N)
title('Neuroticism (N)')
xlabel('Score')
ylabel('# of people')
subplot(4, 2, 4)
histogram(GAD)
title('Generalized Anxiety Disorder (GAD)')
xlabel('Score')
ylabel('# of people')

subplot(4, 2, 5)
scatter(PHQ_prctl, RDS_prctl, 10, 'filled')
xlabel('PHQ')
ylabel('RDS')
title('PHQ vs RDS (Percentiles)')
subplot(4, 2, 6)
scatter(GAD_prctl, N_prctl, 10, 'filled')
xlabel('GAD')
ylabel('N')
title('GAD vs N (Percentiles)')
subplot(4, 2, 7)
scatter(N_prctl, RDS_prctl, 10, 'filled')
xlabel('N')
ylabel('RDS')
title('N vs RDS (Percentiles)')

t = annotation('textbox', [0.09, .93, 0, 0], 'string', 'A')
t.FontWeight = 'bold';
t = annotation('textbox', [0.52, .93, 0, 0], 'string', 'B')
t.FontWeight = 'bold';
t = annotation('textbox', [0.09, .71, 0, 0], 'string', 'C')
t.FontWeight = 'bold';
t=annotation('textbox', [0.52, .71, 0, 0], 'string', 'D')
t.FontWeight = 'bold';
t=annotation('textbox', [0.09, .50, 0, 0], 'string', 'E')
t.FontWeight = 'bold';
t=annotation('textbox', [0.52, .50, 0, 0], 'string', 'F')
t.FontWeight = 'bold';
t=annotation('textbox', [0.09, .29, 0, 0], 'string', 'G')
t.FontWeight = 'bold';

%% grant
figure
sgtitle('Mental Health Questionnaire Comparison')
subplot(1, 2, 1)
scatter(PHQ_prctl, RDS_prctl, 10, 'filled')
xlabel('PHQ')
ylabel('RDS')
title('PHQ vs RDS (Percentiles)')
subplot(1, 2, 2)
scatter(N_prctl, RDS_prctl, 10, 'filled')
xlabel('Neuroticism')
ylabel('RDS')
title('Neuroticism vs RDS (Percentiles)')
t = annotation('textbox', [0.07, .93, 0, 0], 'string', 'A')
t.FontWeight = 'bold';
t = annotation('textbox', [0.5, .93, 0, 0], 'string', 'B')
t.FontWeight = 'bold';
set(gcf,'Position',[471.0000  560.5000  790.5000  338.5000],'PaperPositionMode','auto')
print(gcf,'Fig7','-dpng','-r300');

%%
%get(gcf,'Position')
set(gcf,'Position',[471.0000  560.5000  790.5000  338.5000],'PaperPositionMode','auto')
print(gcf,'Fig7','-dpng','-r300');

%% cronbach for all
N_all = [Mood Miserable Irritable Sensitive Fedup Nervous Worry Tense Embarrassment Nerves Loneliness Guilt];
GAD_all = [Nervous Control Worry Relax Restless Annoyed Afraid];
N_alpha = cronbach(N_all)
GAD_alpha = cronbach(GAD_all)
PHQ_alpha = cronbach(Dep_all)
RDS_alpha = cronbach(D4)

