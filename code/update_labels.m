%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% UPDATE ATLAS LABELS WITH THE REMOVAL OF AREAS WITH TIMECOURSES AT ZERO
% 
% Written by
% Ludovica Romanin - november 2020
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; 

atlas_labels = readtable('../data/atlas/dbs80symm_labels.txt');

load ../data/TCS/areas_zero.mat;

%% 

% remove text labels of areas at zero
atlas_labels(areas_zero,:) = [];

%%

% update the labels' numbers 
n_areas = height(atlas_labels);
new_vars = array2table(1:n_areas);
rows = rows2vars(new_vars);
atlas_labels(:,1) = rows(:,2);

%%

% write 
writetable(atlas_labels,'../data/atlas/dbs80symm_labels_NEW.txt');
