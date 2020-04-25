
load sc_dbs80.mat;
C=sc;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
C(areas_zero,:) = [];
C(:,areas_zero) = [];

load  optimizedhopfawake.mat;

% Optimal G for KL
[M,I] = min(klpstatesControl);
we_optim_kl = WE(I);

% Optimal G for KS
[M1,I1] = min(ksdist);
we_optim_ks = WE(I1);

subplot(1,2,1)
imagesc(log(C));

subplot(1,2,2)
imagesc(log(squeeze(Coptim(I1,:,:))));

figure
hold on;
plot(WE,ksdist,'r');
plot(WE,klpstatesControl,'b');
xline(we_optim_kl,'--b');
xline(we_optim_ks,'--r');
title('Statistics for Control group model')
xlabel('global coupling factor')
legend('KS distance', 'K-L distance')

