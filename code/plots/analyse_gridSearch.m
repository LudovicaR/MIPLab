%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANALYSIS OF GRID SEARCH RESULTS
% 
% Written by
% Ludovica Romanin - november 2020
% 
% This script takes care of plotting the different fitting measures, for a
% and G in the ranges tested in the grid search simulation.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% analysis of the results from the grid search 
% for bifurcation parameter a and global coupling factor G

clear all; clc; 

% load grid search results 
[file,path] = uigetfile(['../../results/gridSearch/','*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/gridSearch/',file]);

%% heatmats showing the different dynamical regimes 

a_coeff = [-0.5:0.02:0.5];

figure()
colormap('jet')
imagesc(WE, a_coeff, fitt')
colorbar
xlabel('Coupling factor (G)')
ylabel('Bifurcation parameter (a)')
title('FC correlation')
%saveas(gcf, '../figures/FC_correlation_symmSC.png')
%%
figure()
colormap('jet')
imagesc(WE, a_coeff, ksdist')
colorbar
xlabel('Coupling factor (G)')
ylabel('Bifurcation parameter (a)')
title('KS distance FCD')
%saveas(gcf, '../figures/KS_distance_FCD_symmSC.png')

%%
figure()
colormap('jet')
imagesc(WE, a_coeff, klpstates')
colorbar
xlabel('Coupling factor (G)')
ylabel('Bifurcation parameter (a)')
title('KL distance PMS')
%saveas(gcf, '../figures/KL_distance_FCD_symmSC.png')

%%
figure()
colormap('jet')
imagesc(WE, a_coeff, entropydist_pms')
colorbar
xlabel('Coupling factor (G)')
ylabel('Bifurcation parameter (a)')
title('Entropy')
%saveas(gcf, '../figures/entropy_symmSC.png')

%%
figure()
colormap('jet')
imagesc(WE, a_coeff, metastability')
colorbar
xlabel('Coupling factor (G)')
ylabel('Bifurcation parameter (a)')
title('Metastability')
%saveas(gcf, '../figures/Metastability_symmSC.png')

%%
idx = find(a_coeff == -0.2);
figure()
plot(WE, fitt(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('FC correlation')
title('a = -0.2')
%saveas(gcf, '../figures/a_02_FCcorr_symmSC.png')

%%
idx = find(a_coeff == 0.0);
figure()
plot(WE, fitt(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('FC correlation')
title('a = 0.0')
%saveas(gcf, '../figures/a00_FCcorr_symmSC.png')

%%
idx = find(a_coeff == 0.2);
figure()
plot(WE, fitt(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('FC correlation')
title('a = 0.2')
%saveas(gcf, '../figures/a02_FCcorr_symmSC.png')

%%
idx = find(a_coeff == -0.2);
figure()
plot(WE, ksdist(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('KS distance FCD')
title('a = -0.2')
%saveas(gcf, '../figures/a_02_KS_symmSC.png')

%%
idx = find(a_coeff == 0.0);
figure()
plot(WE, ksdist(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('KS distance FCD')
title('a = 0.0')
%saveas(gcf, '../figures/a00_KS_symmSC.png')

%%
idx = find(a_coeff == 0.2);
figure()
plot(WE, ksdist(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('KS distance FCD')
title('a = 0.2')
%saveas(gcf, '../figures/a02_KS_symmSC.png')

%%
idx = find(a_coeff == -0.2);
figure()
plot(WE, metastability(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('Metastability')
title('a = -0.2')
%saveas(gcf, '../figures/a_02_meta_symmSC.png')

%%
idx = find(a_coeff == 0.0);
figure()
plot(WE, metastability(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('Metastability')
title('a = 0.0')
%saveas(gcf, '../figures/a00_meta_symmSC.png')

%%
idx = find(a_coeff == 0.2);
figure()
plot(WE, metastability(idx,:), '-b')
xlabel('Coupling factor (G)')
ylabel('Metastability')
title('a = 0.2')
%saveas(gcf, '../figures/a02_meta_symmSC.png')

%%
best_kl = min(klpstates(:));
best_ks = min(ksdist(:)); 

[idx_best_kl_1, idx_best_kl_2] = find(klpstates == best_kl);
[idx_best_ks_1, idx_best_ks_2] = find(ksdist == best_ks);

%%
figure()
hold on
plot(WE, fitt(idx_best_kl_1,:), '-k')
plot(WE, fitt(idx_best_ks_1,:), '-c')
xline(WE(idx_best_kl_2),'--b');
xline(WE(idx_best_ks_2),'--r');
legend('Pearson corr.  - KL', 'Pearson corr.  - KS', 'KL optim', 'KS optim')
xlabel('coupling factor G')

figure()
hold on
plot(a_coeff, fitt(:,idx_best_kl_2), '-k')
plot(a_coeff, fitt(:,idx_best_ks_2), '-c')
xline(a_coeff(idx_best_kl_1),'--b');
xline(a_coeff(idx_best_ks_1),'--r');
legend('Pearson corr.  - KL', 'Pearson corr.  - KS', 'KL optim', 'KS optim')
xlabel('Bifurcation parameter (a)')

%%
FCsim_optim_kl = squeeze(FCsim(idx_best_kl_1,idx_best_kl_2,:,:));
FCsim_optim_ks = squeeze(FCsim(idx_best_ks_1, idx_best_ks_2,:,:));

%%
figure()
colormap(jet)
hold on
subplot(1,3,1)
imagesc(FCsim_optim_kl)
colorbar
title('FC sim - KL')

subplot(1,3,2)
imagesc(FCemp)
colorbar
title('FC emp')

subplot(1,3,3)
imagesc(FCsim_optim_ks)
colorbar
title('FC sim - KS')
