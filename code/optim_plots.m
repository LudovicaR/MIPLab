
clear all;
%load sc_dbs80.mat;
%C=sc;

load meanSC_56HC_Desikan.mat
SC_56 = meanSC;
load meanSC_56HC_Desikan_woCC.mat
SC_woCC = meanSC;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
SC_56(areas_zero,:) = [];
SC_56(:,areas_zero) = [];

load  optimizedhopfawake_56HC.mat;
load empiricalLEiDA.mat;


%% find optimal G and corresponding PMS

% Optimal G for KL
[M,I] = min(klpstates);
we_optim_kl = WE(I);

% Optimal G for KS
[M1,I1] = min(ksdist);
we_optim_ks = WE(I1);

figure
hold on;
plot(WE,ksdist,'r');
plot(WE,klpstates,'b');
xline(we_optim_kl,'--b');
xline(we_optim_ks,'--r');
title('Statistics for Control model')
xlabel('global coupling factor')
legend('KS distance', 'KL distance')

figure
plot(WE, fitt, 'k');
xline(we_optim_kl,'--b');
xline(we_optim_ks,'--r');
xlabel('global coupling factor')
legend('Pearson correlation', 'G_{opt}, KL', 'G_{opt}, KS', 'FontSize', 11)

%% SC matrices plot
figure
subplot(1,2,1)
imagesc(log(SC_56)); 
colorbar
title('healthy SC')
hold on 

subplot(1,2,2)
imagesc(log(SC_woCC));
colorbar
title('lesioned SC (w/o CC)')
hold on

