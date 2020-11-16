%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANALYSIS OF THE MODEL PARAMETERS BASED ON EMPIRICAL DATA
%
% Written by
% Ludovica Romanin - november 2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; 

%load data matrices for the two groups
load AgCC_TCS_new.mat
load Controls_TCS.mat

%%
%parameters
N_areas=size(Controls_TCS, 1);
TP=size(Controls_TCS, 2);
NSUB_Controls=size(Controls_TCS,3);
NSUB_AgCC=size(AgCC_TCS,3);

%% parameters 
TR=2;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

%% FC matric - Control 

R_subjs_ctrl = zeros(NSUB_Controls, TP, TP);
P_subjs_ctrl = zeros(NSUB_Controls, TP, TP);
mean_fc_subjs_ctrl = zeros(NSUB_Controls, TP);
std_fc_subjs_ctrl = zeros(NSUB_Controls, TP);
k_subjs_ctrl = zeros(NSUB_Controls, TP);

for s=1:NSUB_Controls
    BOLD = Controls_TCS(:,:,s);
    [R,P] = corrcoef(BOLD);
    R_subjs_ctrl(s,:,:) = R;
    P_subjs_ctrl(s,:,:) = P;
    mean_fc_subjs_ctrl(s,:) = mean(R);
    std_fc_subjs_ctrl(s,:) = std(R);

    Phase_BOLD=zeros(N_areas,TP);
    
    % Get the BOLD phase using the Hilbert transform
    for seed=1:N_areas
        BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
        signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end
    
    % Kuramoto order parameter --> measure of phase locking 
    k = zeros(1,TP);
    for t = 1:TP
        phase_sum = 0;
        for seed=1:N_areas
            phase_sum = phase_sum + cos(Phase_BOLD(seed,t))+1i*sin(Phase_BOLD(seed,t));
        end
        k(t) = abs(1/N_areas *phase_sum);
    end
    k_subjs_ctrl(s,:) = k;
end

%% FC matric - AgCC 

R_subjs_agcc = zeros(NSUB_AgCC, TP, TP);
P_subjs_agcc = zeros(NSUB_AgCC, TP, TP);
mean_fc_subjs_agcc = zeros(NSUB_AgCC, TP);
std_fc_subjs_agcc = zeros(NSUB_AgCC, TP);
k_subjs_agcc = zeros(NSUB_AgCC, TP);

for s=1:NSUB_AgCC
    BOLD = AgCC_TCS(:,:,s);
    [R,P] = corrcoef(BOLD);
    R_subjs_agcc(s,:,:) = R;
    P_subjs_agcc(s,:,:) = P;
    mean_fc_subjs_agcc(s,:) = mean(R);
    std_fc_subjs_agcc(s,:) = std(R);
    
    Phase_BOLD=zeros(N_areas,TP);
    
    % Get the BOLD phase using the Hilbert transform
    for seed=1:N_areas
        BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
        signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end
    
    % Kuramoto order parameter --> measure of phase locking 
    k = zeros(1,TP);
    for t = 1:TP
        phase_sum = 0;
        for seed=1:N_areas
            phase_sum = phase_sum + cos(Phase_BOLD(seed,t))+1i*sin(Phase_BOLD(seed,t));
        end
        k(t) = abs(1/N_areas *phase_sum);
    end
    k_subjs_agcc(s,:) = k;
    
end

%% plot of level of synchronization through Kuramoto parameter 

figure()
hold on
stdshade(k_subjs_ctrl,0.4,'b')
stdshade(k_subjs_agcc,0.4,'r')
legend('std ctrl', 'mean ctrl', 'std agcc', 'mean agcc')

%%
% analysis of fitt = pearson correlation 
% analysis of metastability = std of phase synchrony 
clear all; clc; 
load("../results/data_FCD/HopfModel_results_56HC.mat");
fitt_ctrl = fitt;
C_ctrl = sumC; 
metastability_ctrl = metastability;

load("../results/data_FCD/HopfModel_results_56HC_woCC.mat");
fitt_wocc = fitt;
C_wocc = sumC; 
metastability_wocc = metastability;

figure()
hold on
plot(WE, fitt_ctrl, '-r')
plot(WE, fitt_wocc, '-b')
legend('ctrl', 'wocc')

C_diff = C_ctrl-C_wocc;
C_diff(:,1)

figure()
hold on
plot(WE, metastability_ctrl, '-r')
plot(WE, metastability_wocc, '-b')
plot(WE, ones(1,NWE)*metastabilitydata, '-k')
legend('ctrl', 'wocc')


%%
clear all; clc; 
load AgCC_TCS_new.mat
load Controls_TCS.mat
load("../results/empiricalLEiDA_K_10_15.mat");
load("../results/LEiDA_results_K_10_15.mat");

%parameters
N_areas=size(Controls_TCS, 1);
TP=size(Controls_TCS, 2);
NSUB_Controls=size(Controls_TCS,3);
NSUB_AgCC=size(AgCC_TCS,3);

c_agcc = zeros(NSUB_AgCC,TP,TP);
c_ctrl = zerso(NSUB_Controls,TP,TP);

for s=1:NSUB_AgCC
    c_agcc(s,:,:) = cov(AgCC_TCS(:,:,s));
    
end

for s=1:NSUB_Controls
    c_ctrl(s,:,:) = cov(Controls_TCS(:,:,s));
end

%% synchrony per subject 
figure
hold on
subplot(1,2,1)
bar(mean(k_subjs_ctrl,2),'b', 'FaceAlpha',0.2)
hold on
errorbar(mean(k_subjs_ctrl,2),std(k_subjs_ctrl,0,2), 'Color', 'b', 'LineStyle', 'none')
legend('ctrl')
subplot(1,2,2)
bar(mean(k_subjs_agcc,2),'r', 'FaceAlpha',0.2)
hold on 
errorbar(mean(k_subjs_agcc,2),std(k_subjs_agcc,0,2), 'Color', 'r','LineStyle', 'none')
legend('agcc')
