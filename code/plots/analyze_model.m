%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WHOLE BRAIN CONNECTIVITY MEASURES
%
% Written by
% Ludovica Romanin - november 2020
% 
% This script calculates the different connectivity measures (based on the
% Kuramoto Order Parameter (KOP), to investigate the level of
% synchronization.
% For more information refer to these publications: 
% (Deco, 2017) The dynamics of resting fluctuations in the brain 
%   (doi:10.1038/s41598-017-03073-5)
% 
% (Demirtas, 2017) whole-brain modeling of Alzheimer's disease 
%   (doi:10.1016/j.nicl.2017.08.006)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; 

% add functions folder to the path
addpath('../functions/')

%load data matrices for the two groups
load ../../data/TCS/AgCC_TCS.mat
load ../../data/TCS/Controls_TCS.mat


%%
%parameters
TP=size(Controls_TCS, 2); % number of timepoints 
NSUB_Controls=size(Controls_TCS,3);  % NSUB = number of subjects 
NSUB_AgCC=size(AgCC_TCS,3);

%% remove areas at zero to avoid NaN
load ../../data/TCS/areas_zero.mat;
Controls_TCS(areas_zero,:,:) = [];
AgCC_TCS(areas_zero,:,:) = [];

N_areas=size(Controls_TCS, 1); % final number of areas 

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

%% FC computation = Pearson'correlation between timeseries of brain regions

% Control subjects
FC_ctrl = zeros(NSUB_Controls, N_areas, N_areas); 
for s=1:NSUB_Controls
    BOLD = Controls_TCS(:,:,s);
    FC_ctrl(s,:,:) = corrcoef(BOLD');
end

meanFC_ctrl = squeeze(mean(FC_ctrl));

% AgCC subjects
FC_agcc = zeros(NSUB_AgCC, N_areas, N_areas);
for s=1:NSUB_AgCC
    BOLD = AgCC_TCS(:,:,s);
    FC_agcc(s,:,:) = corrcoef(BOLD');
end

meanFC_agcc = squeeze(mean(FC_agcc));

%% Pearson correlation between FC of the two groups (agcc + ctrl)
Isubdiag = find(tril(ones(N_areas),-1));
cc = corrcoef(atanh(meanFC_ctrl(Isubdiag)),atanh(meanFC_agcc(Isubdiag)));
FC_corr = cc(2);

clear cc Isubdiag

%% calculate FC strength, measure of regional connectivity 
FCstrength_ctrl = sum(meanFC_ctrl);
FCstrength_agcc = sum(meanFC_agcc);

%% Kuramoto order paremeter (KOP) - Control
K_ctrl = zeros(NSUB_Controls,TP);
for s=1:NSUB_Controls
    % Get the BOLD signals from this subject in this task
    BOLD = Controls_TCS(:,:,s);
    Phase_BOLD=zeros(N_areas,TP);
    
    % Get the BOLD phase using the Hilbert transform
    for seed=1:N_areas
        BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
        signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end
    
    for t=1:TP
        
        %Calculate the Instantaneous FC (BOLD Phase Synchrony)
        iFC=zeros(N_areas);
        for n=1:N_areas
            for p=1:N_areas
                iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end
        
        phase_sum = 0;
        for seed=1:N_areas
            % Euler formula:exp(jt) = cos(t) + jsin(t)
            phase_sum = phase_sum + cos(iFC(seed))+1i*sin(iFC(seed));
        end
        
        K_ctrl(s,t) = 1/N_areas * phase_sum;
        
    end
end

%% Kuramoto order paremeter (KOP) - AgCC
K_agcc = zeros(NSUB_AgCC,TP);
for s=1:NSUB_AgCC
    % Get the BOLD signals from this subject in this task
    BOLD = AgCC_TCS(:,:,s);
    Phase_BOLD=zeros(N_areas,TP);
    
    % Get the BOLD phase using the Hilbert transform
    for seed=1:N_areas
        BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
        signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
    end
    
    for t=1:TP
        
        %Calculate the Instantaneous FC (BOLD Phase Synchrony)
        iFC=zeros(N_areas);
        for n=1:N_areas
            for p=1:N_areas
                iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end
        
        phase_sum = 0;
        for seed=1:N_areas
            % Euler formula:exp(jt) = cos(t) + jsin(t)
            phase_sum = phase_sum + cos(iFC(seed))+1i*sin(iFC(seed));
        end
        
        K_agcc(s,t) = 1/N_areas * phase_sum;
        
    end
end

%% Global coherence and metastability 

% global coherence or mean synchronization = temporal average of KOP
sync_ctrl = mean(K_ctrl,2);
sync_agcc = mean(K_agcc,2);

% metastability (variation in synchronization over time) = standard deviation of KOP
metastability_ctrl = std(K_ctrl,[],2);
metastability_agcc = std(K_agcc,[],2);

%% plot of level of synchronization through Kuramoto parameter 

figure()
hold on
stdshade(K_agcc,0.4,'b')
%stdshade(k_subjs_agcc,0.4,'r')
xlabel('Time (seconds)','FontSize', 14)
ylabel('KOP', 'FontSize', 14)
%legend('std ctrl', 'mean ctrl', 'std agcc', 'mean agcc')


%% average synchrony across subjects
figure
hold on
bar(1, median(sync_ctrl),'b', 'FaceAlpha',0.2)
hold on
errorbar(1, median(sync_ctrl),std(sync_ctrl), 'Color', 'b', 'LineStyle', 'none')

bar(2, median(sync_agcc),'r', 'FaceAlpha',0.2)
hold on 
errorbar(2, median(sync_agcc),std(sync_ctrl), 'Color', 'r','LineStyle', 'none')
xticks([1, 2])
xticklabels({'Control', 'AgCC'})
title('global coherence - mean synchronization')

