clear all; close all; clc

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/')

% Run CCA on all subjects
fprintf('Running CCA on all subjects \n');
load('CCA_inputs_behavior_CCA.mat');
load('CCA_inputs_brain_CCA.mat');

%[pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],Dall,2000);
%save('CCA_results_Dall.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V
[pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT_NEW],Dsum,2000);
save('CCA_results_DsumNEW.mat','pfwer','r','A','B','U','V'); 

% Run CCA on patients only
% fprintf('Running CCA on patients only \n');
% load('CCA_inputs_behavior_CCApatients.mat');
% load('CCA_inputs_brain_CCApatients.mat');
%
% [pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],Dall,2000);
% save('CCA_results_patients_Dall.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V
% [pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],Dsum,2000);
% save('CCA_results_patients_Dsum.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V

% Run CCA for late life depression only
% fprintf('Running CCA on late life depression only \n');
% load('CCA_inputs_behavior_CCAlatelife.mat');
% load('CCA_inputs_brain_CCAlatelife.mat');
% 
% %[pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],Dall,2000);
% %save('CCA_results_latelife_Dall.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V
% %[pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],Dsum,2000);
% %save('CCA_results_latelife_Dsum.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V
% [pfwer,r,A,B,U,V] = permcca([MAPS_OUT REST_OUT IDP_OUT],D4,2000);
% save('CCA_results_latelife_D4.mat','pfwer','r','A','B','U','V'); clear pfwer r A B U V
