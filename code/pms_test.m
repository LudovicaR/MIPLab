clear all;
load meanSC_56HC_Desikan.mat;
C=meanSC;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
C(areas_zero,:) = [];
C(:,areas_zero) = [];

load  optimizedhopfawake_56HC.mat;
load empiricalLEiDA.mat P1emp P2emp PTR1emp PTR2emp;

% Optimal G for KL
[M,I] = min(klpstates);
we_optim_kl = WE(I);

% Optimal G for KS
[M1,I1] = min(ksdist);
we_optim_ks = WE(I1);

%% plot PMS, simulated and empirical

% optimizing on KS 
pms = Pstatessimul(I1,:);
figure
subplot(1,3,1)
bar(pms)
title('simulated (KS)')
subplot(1,3,2)
mean_pms_control = mean(P2emp,1);
bar(mean_pms_control)
title('control')
subplot(1,3,3)
mean_pms_agCC = mean(P1emp,1);
bar(mean_pms_agCC)
title('AgCC')


% optimizing on KL
pms = Pstatessimul(I,:);
figure
subplot(1,3,1)
bar(pms)
title('simulated (KL)')
subplot(1,3,2)
mean_pms_emp = mean(P2emp,1);
bar(mean_pms_emp)
title('control')
subplot(1,3,3)
mean_pms_agCC = mean(P1emp,1);
bar(mean_pms_agCC)
title('AgCC')


%% t_test to assess significance of difference between PMS

[H_control,P_control, CI_control, STATS_control] = ttest2(Pstatessimul(I,:).*ones(28,1), P2emp);
[H_agcc, P_agcc, CI_agcc, STATS_agcc] = ttest2(Pstatessimul(I,:).*ones(13,1),P1emp);

%%
figure
subplot(1,3,1)
imagesc(PTR1emp)
colorbar
subplot(1,3,2)
imagesc(PTR2emp)
colorbar
subplot(1,3,3)
imagesc(squeeze(PTRsimul(I,:)));
colorbar
