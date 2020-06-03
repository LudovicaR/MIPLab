
clear all;
load meanSC_56HC_Desikan_woCC.mat;
C=meanSC;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
C(areas_zero,:) = [];
C(:,areas_zero) = [];

load  optimizedhopfawake_56HC_woCC.mat;
load empiricalLEiDA.mat;

load HopfModel_results_56HC_woCC.mat FC_simul FCemp FCphasesemp;

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
title('Statistics for Virtual Callosotomy')
xlabel('global coupling factor')
legend('KS distance', 'KL distance')

figure
plot(WE, fitt, 'k');
xline(we_optim_kl,'--b');
xline(we_optim_ks,'--r');
xlabel('global coupling factor')
legend('Pearson correlation', 'G_{opt}, KL', 'G_{opt}, KS', 'FontSize', 11)

%% plot FC matrices
figure
subplot(1,2,1)
imagesc(FC_simul)
title('simulated')
hold on

subplot(1,2,2)
imagesc(FCemp)
title('empirical')
hold on

%% FCD plots 
load fcd_result_56HC_woCC.mat

%%
fcd_kl = phfcd_sim(I,:,:);
fcd_ks = phfcd_sim(I1,:,:);

fcd_kl = squeeze(fcd_kl);
fcd_ks = squeeze(fcd_ks);

%%
N = 74;
TP = 200;
pat_size = N*(N-1) + (TP-21) - 2;
fcd_length = pat_size*(pat_size-1)/2;

fcd_length_area = floor(fcd_length/74);

%%

fcd_kl_areas = zeros(N,fcd_length_area);
fcd_ks_areas = zeros(N,fcd_length_area);

for seed=1:N
    fcd_kl_areas(seed,:) = fcd_kl((seed-1)*fcd_length_area+1 : seed*fcd_length_area);
    fcd_ks_areas(seed,:) = fcd_ks((seed-1)*fcd_length_area+1 : seed*fcd_length_area);
end

%%
std_fcd_kl = std(fcd_kl_areas,0,2);
std_fcd_ks = std(fcd_ks_areas,0,2);

std_kl = std(fcd_kl);
std_ks = std(fcd_ks);
std_emp = std(phfcddata_emp);

%%
atlas_labels = readtable('dbs80symm_labels_NEW.txt');
aal =  atlas_labels.Var2;
%%
x = linspace(1,74,74);
std_bar = [std_fcd_kl, std_fcd_ks];
b = bar(x,std_bar);