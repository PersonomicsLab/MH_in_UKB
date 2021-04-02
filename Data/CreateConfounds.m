clear all; close all; clc

% Based on https://git.fmrib.ox.ac.uk/falmagro/ukb_unconfound_v2/-/blob/master/conf_processing/common_matlab/duplicateDemedianNormBySite.m
% and https://git.fmrib.ox.ac.uk/falmagro/ukb_unconfound_v2/-/blob/master/conf_processing/common_matlab/duplicateCategorical.m

addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts/FSLNets')
addpath('/scratch/janine/MentalHealthInUKB/MatlabScripts')
INPUT = 'Subjects_CCA.csv';

%% Load confound information
S = load(sprintf('%s/%s','/scratch/janine/MentalHealthInUKB/SubjectSplits/',INPUT));
Conf = readtable('Conf_scan1.tsv','FileType','text');
[~,s,~] = intersect(table2array(Conf(:,1)),S); Conf = Conf(s,:);
H = get_UKB_headers(Conf);

%% Create site varriables
finalConfs = [];
n1 = strfind(H,'54-2.0'); n1 = find(~cellfun(@isempty,n1));
site = table2array(Conf(:,n1)); clear n1
values = site;
diffValues = unique(values);
numDiffValues = length(diffValues);
if numDiffValues > 1
    numNewCols = numDiffValues-1;
    newConfound = zeros(length(S), numDiffValues-1);
    
    % Subjects (of this site) with the first value get a -1 in all new columns
    indFinal = find(values == diffValues(1));
    newConfound(indFinal,:)=-1;
    
    % Each new value gets a new column. Subjects (of this Site) with
    % this value have 1 in this column. All other subjects have a 0.
    for j=2:numDiffValues
        indFinal = find(values == diffValues(j));
        newConfound(indFinal,j-1)=1;
        
        indNotZero = find(newConfound(:,j-1) ~= 0);
        tmpVar = newConfound(indNotZero,j-1);
        
        % Normalising only makes sense if there is more than 1
        % different value
        newConfound(indNotZero,j-1) = nets_normalise(newConfound(indNotZero,j-1));
    end
    
end

%% Creating other variables
% Age
n1 = strfind(H,'21003-2.0'); n1 = find(~cellfun(@isempty,n1));
values = table2array(Conf(:,n1)); clear n1
% Age squared
values(:,2) = values(:,1).^2;
% Sex
n1 = strfind(H,'31-0.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Age * sex
values(:,4) = values(:,1).*values(:,3);
% Head size
n1 = strfind(H,'25000-2.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Head motion rfMRI
n1 = strfind(H,'25741-2.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Head motion tfMRI
n1 = strfind(H,'25742-2.0'); n1 = find(~cellfun(@isempty,n1));
values = [values table2array(Conf(:,n1))]; clear n1
% Date
n1 = strfind(H,'53-2.0'); n1 = find(~cellfun(@isempty,n1));
days = table2array(Conf(:,n1)); clear n1
days = datenum(days); days = days - (min(days)-1);
values = [values days]; clear days
% Date squared
values(:,9) = values(:,8).^2;

%% Split regressors for sites
numVars  = size(values,2);
subjectIndicesBySite{1} = find(site==11025);
subjectIndicesBySite{2} = find(site==11027);
numSites = length(subjectIndicesBySite);

% Demedian globally each column
for i = 1 : numVars
    values(:,i) = values(:,i) - nanmedian(values(:,i));
    madV = mad(values(:,i),1) * 1.48;
    if madV < eps
        madV = nanstd(values(:,i));
    end
    values(:,i) = values(:,i) / madV;
end

% Split the colums: 1 per site. 0-padding by default.
finalConfs = [finalConfs zeros(length(S),numSites*numVars)];

% For each Site
for i = 1:numSites
    % For each variable (column) that is received in "values"
    for j = 1:numVars
        V = values(subjectIndicesBySite{i},j);
        
        medTP = nanmedian(V);
        
        V(find(V > 8)) =  NaN;
        V(find(V <-8)) =  NaN;
        
        V(isnan(V)) = medTP;
        
        V = nets_normalise(V);
        finalConfs(subjectIndicesBySite{i}, (numVars*(i-1) + j)) = V;
    end
end

conf=[finalConfs, newConfound];
save(sprintf('Confounds_Subjects_%s.mat',INPUT(10:end-4)),'conf');
dlmwrite(sprintf('Confounds_Subjects_%s.csv',INPUT(10:end-4)),conf,'precision',25);
