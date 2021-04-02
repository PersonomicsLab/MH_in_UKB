clear all; close all; clc

%%% RUN ON MY OWN LAPTOP %%%

% Load phenotype variable IDs
load vars.mat

% Write confound text files for scan 1
system('rm Conf_scan1.txt');
Idone = [];
fileID = fopen('Conf_scan1.txt','w');
I = find(contains(Conf(:,2),'Sex'));
fprintf(fileID,'%s-0.0\n',Conf(I,1)); Idone = [Idone I]; clear I
I = find(contains(Conf(:,2),'MH'));
fprintf(fileID,'%s-0.0\n',Conf(I,1)); Idone = [Idone I]; clear I
for n = setdiff(1:length(Conf),Idone)
    fprintf(fileID,'%s-2.0\n',Conf(n,1));
end
fclose(fileID); clear I n Idone fileID

% Write confound text files for scan 2
system('rm Conf_scan2.txt');
Idone = [];
fileID = fopen('Conf_scan2.txt','w');
I = find(contains(Conf(:,2),'Sex'));
fprintf(fileID,'%s-0.0\n',Conf(I,1)); Idone = [Idone I]; clear I
I = find(contains(Conf(:,2),'MH'));
fprintf(fileID,'%s-0.0\n',Conf(I,1)); Idone = [Idone I]; clear I
for n = setdiff(1:length(Conf),Idone)
    fprintf(fileID,'%s-3.0\n',Conf(n,1));
end
fclose(fileID); clear I n Idone fileID

% Write MH text file for scan 1
system('rm MH_scan1.txt');
fileID = fopen('MH_scan1.txt','w');
for n = 1:length(MH)
    if strlength(MH{n,1})==4
        fprintf(fileID,'%s-2.0\n',MH(n,1));
    elseif strlength(MH{n,1})==5
        fprintf(fileID,'%s-0.0\n',MH(n,1));
    end
end
fclose(fileID); clear n fileID

% Write MH text file for scan 2
system('rm MH_scan2.txt');
fileID = fopen('MH_scan2.txt','w');
for n = 1:length(MH)
    if strlength(MH{n,1})==4
        fprintf(fileID,'%s-3.0\n',MH(n,1));
    elseif strlength(MH{n,1})==5
        fprintf(fileID,'%s-0.0\n',MH(n,1));
    end
end
fclose(fileID); clear n fileID

% Write IDP text file for scan 1
system('rm IDP_scan1.txt');
fileID = fopen('IDP_scan1.txt','w');
for n = 1:length(IDP)
    fprintf(fileID,'%s-2.0\n',IDP(n,1));
end
fclose(fileID); clear n fileID

% Write IDP text file for scan 1
system('rm IDP_scan2.txt');
fileID = fopen('IDP_scan2.txt','w');
for n = 1:length(IDP)
    fprintf(fileID,'%s-3.0\n',IDP(n,1));
end
fclose(fileID); clear n fileID

% Write extras file
system('rm Extras.txt');
fileID = fopen('Extras.txt','w');
for n = 1:length(IDP)
    fprintf(fileID,'20433-0.0\n'); % Age at first episode of depression
end
fclose(fileID); clear n fileID

% Run phenotype extraction
%system('funpack -q -ow -s subs20k.csv -co Conf_scan1.txt Conf_scan1.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
%system('funpack -q -ow -s subs20k.csv -co Conf_scan2.txt Conf_scan2.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
%system('funpack -q -ow -s subs20k.csv -co MH_scan1.txt MH_scan1.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
%system('funpack -q -ow -s subs20k.csv -co MH_scan2.txt MH_scan2.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
%system('funpack -q -ow -s subs20k.csv -co IDP_scan1.txt IDP_scan1.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
%system('funpack -q -ow -s subs20k.csv -co IDP_scan2.txt IDP_scan2.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
system('funpack -q -ow -s subs20k.csv -co Extras.txt Extras.tsv /Users/janinebijsterbosch/Dropbox/WashU/Data/UKB/pheno/ukb41049.csv');
