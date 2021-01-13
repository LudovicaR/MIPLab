%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create structure with all the timecourses
%
% Written by
% Ludovica Romanin - november 2020
% 
% Final output dimension: Nsubjs x Nareas x Timepoints
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

tic
% AgCC

path = '../data/246BN/AgCC';

files = dir(fullfile(path, '*.mat'));
tc_test_size = load(fullfile(files(1).folder, files(1).name));
AgCC_TCS = zeros(size(tc_test_size.TCS,1), size(tc_test_size.TCS,2),length(files));
for i=1:length(files)
    tc = load(fullfile(files(i).folder, files(i).name));
    tc = tc.TCS;
    AgCC_TCS(:,:,i) = tc;
end

%save('../data/TCS/AgCC_TCS_.mat', 'AgCC_TCS'); 

% CONTROL

path = '../data/246BN/HC';

files = dir(fullfile(path, '*.mat'));
tc_test_size = load(fullfile(files(1).folder, files(1).name));
Controls_TCS = zeros(size(tc_test_size.TCS,1), size(tc_test_size.TCS,2),length(files));
for i=1:length(files)
    tc = load(fullfile(files(i).folder, files(i).name));
    tc = tc.TCS;
    Controls_TCS(:,:,i) = tc;
end

%save('../data/TCS/Controls_TCS_.mat', 'Controls_TCS'); 

elapsed_time = toc;
toc