clear all; close all; clc

addpath('/Users/janinebijsterbosch/Box/00_CHPC2_Backups/MentalHealthInUKB_CHPC/MatlabScripts/')

% Run CCA on all subjects
fprintf('Running CCA on all subjects \n');
load('CCA_inputs_behavior_CCAnew.mat');
load('CCA_inputs_brain_CCAnew.mat');

[pfwer,r,A,B,U,V] = permcca([REST_OUT STRUCT_OUT TASK_OUT],Dsum,2000);
save('CCA_results_CCAnew.mat','pfwer','r','A','B','U','V'); 

