
% script responsible to create the structure with all the timecourses 
% the final output is of dimension: Nsubjs x Nareas x Timepoints

clear all; clc;

tic
% AgCC

path = '/Users/ludovica/Documents/GitHub/MIPLab/data/800Schaefer/AgCC';

files = dir(fullfile(path, '*.mat'));
tc_test_size = load(fullfile(files(1).folder, files(1).name));
AgCC_TCS = zeros(size(tc_test_size.TCSnf,1), size(tc_test_size.TCSnf,2),length(files));
for i=1:length(files)
    tc = load(fullfile(files(i).folder, files(i).name));
    tc = tc.TCSnf;
    AgCC_TCS(:,:,i) = tc;
end

save('AgCC_TCS_800.mat', 'AgCC_TCS'); 

% CONTROL

path = '/Users/ludovica/Documents/GitHub/MIPLab/data/800Schaefer/HC';

files = dir(fullfile(path, '*.mat'));
tc_test_size = load(fullfile(files(1).folder, files(1).name));
Controls_TCS = zeros(size(tc_test_size.TCSnf,1), size(tc_test_size.TCSnf,2),length(files));
for i=1:length(files)
    tc = load(fullfile(files(i).folder, files(i).name));
    tc = tc.TCSnf;
    Controls_TCS(:,:,i) = tc;
end

save('Controls_TCS_800.mat', 'Controls_TCS'); 

elapsed_time = toc;
toc