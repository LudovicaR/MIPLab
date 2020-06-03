atlas_labels = readtable('dbs80symm_labels.txt');

load areas_zero.mat;

%% 

atlas_labels(areas_zero,:) = [];

%%
n_areas = height(atlas_labels);
new_vars = array2table(1:n_areas);
rows = rows2vars(new_vars);
atlas_labels(:,1) = rows(:,2);
%%
writetable(atlas_labels,'dbs80symm_labels_NEW.txt');