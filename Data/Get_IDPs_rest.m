clear alll; close all; clc

subs = load('../SubjectSplits/Subjects_testretest.csv');
Amp100 = zeros(length(subs),55);
Fnet100 = zeros(length(subs),1485);
Pnet100 = zeros(length(subs),1485);
Amp25 = zeros(length(subs),21);
Fnet25 = zeros(length(subs),210);
Pnet25 = zeros(length(subs),210);
Dir = '/NRG-data/NRG/beegfs/mirrir/biobank/derivatives/rfMRI_IDPs/';

for s = 1:length(subs)
    if exist(sprintf('%s/AMP100/%d_25755_3_0.txt',Dir,subs(s)),'file')
    Amp100(s,:) = load(sprintf('%s/AMP100/%d_25755_3_0.txt',Dir,subs(s)));
    end
    if exist(sprintf('%s/FNET100/%d_25751_3_0.txt',Dir,subs(s)),'file')
    Fnet100(s,:) = load(sprintf('%s/FNET100/%d_25751_3_0.txt',Dir,subs(s)));
    end
    if exist(sprintf('%s/PNET100/%d_25753_3_0.txt',Dir,subs(s)),'file')
    Pnet100(s,:) = load(sprintf('%s/PNET100/%d_25753_3_0.txt',Dir,subs(s)));
    end
    if exist(sprintf('%s/AMP25/%d_25754_3_0.txt',Dir,subs(s)),'file')
    Amp25(s,:) = load(sprintf('%s/AMP25/%d_25754_3_0.txt',Dir,subs(s)));
    end
    if exist(sprintf('%s/FNET25/%d_25750_3_0.txt',Dir,subs(s)),'file')
    Fnet25(s,:) = load(sprintf('%s/FNET25/%d_25750_3_0.txt',Dir,subs(s)));
    end
    if exist(sprintf('%s/PNET25/%d_25752_3_0.txt',Dir,subs(s)),'file')
    Pnet25(s,:) = load(sprintf('%s/PNET25/%d_25752_3_0.txt',Dir,subs(s)));
    end
end

IDP = [Amp100 Fnet100 Pnet100 Amp25 Fnet25 Pnet25];
save('IDPs_scan2.mat','IDP','subs');