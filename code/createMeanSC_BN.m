%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create mean structural connectome from all subject's SCs
%
% Written by
% Ludovica Romanin - november 2020
% 
% Script used to obtain the meanSC.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;  

%% load files with structural connectomes - Controls 

path = '../data/SCs_HC_norm';

files = dir(fullfile(path, '*.mat'));
sc_size = load(fullfile(files(1).folder, files(1).name));
meanSC_all = zeros(size(sc_size.Cnorm,1), size(sc_size.Cnorm,2),length(files));
for i=1:length(files)
    sc = load(fullfile(files(i).folder, files(i).name));
    sc = sc.Cnorm;
    meanSC_all(:,:,i) = sc;
end

meanSC = mean(meanSC_all,3);

%save('../data/meanSCs/meanSC_BN_HC.mat', 'meanSC'); 

%% load files with structural connectomes - AgCC 

clear all; clc; 

path = '../data/SCs_AgCC_norm/complete';

subjects = {'s020', 's021', 's022', 's102', 's103', 's105', 's106', 's113'};

files = dir(fullfile(path, '*.mat'));
sc_size = load(fullfile(files(1).folder, files(1).name));
meanSC_all = zeros(size(sc_size.Cnorm,1), size(sc_size.Cnorm,2),length(files));
for i=1:length(files)
    sc = load(fullfile(files(i).folder, files(i).name));
    sc = sc.Cnorm;
    meanSC_all(:,:,i) = sc;
end

meanSC = mean(meanSC_all,3);

%save('../data/meanSCs/meanSC_BN_AgCC.mat', 'meanSC'); 
