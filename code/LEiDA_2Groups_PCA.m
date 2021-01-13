%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS - use PCA instead of k-means
%
% Written by
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
% Giulia Preti March 2020 - adaptation for AgCC analysis: 2 groups
% Ludovica Romanin March 2020 - remove small areas giving timecourses=0
% Ludovica Romanin December 2020 - apply PCA instead of k-means (clusters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

addpath('functions/')

%% Load data matrices for the two groups
% LOAD Controls DATA

% Controls_TCS.mat : parcellation with Desikan atlas, 80 regions 
% Controls_TCS_246BN.mat : parcellation with Brainnetome atlas, 246 regions

load('../data/TCS/Controls_TCS.mat');

% LOAD AgCC DATA

% AgCC_TCS.mat : parcellation with Desikan atlas, 80 regions
% AgCC_TCS_246BN.mat : parcellation with Brainnetome atlas, 246 regions

load('../data/TCS/AgCC_TCS.mat');

%%

%parameters
N_areas=size(Controls_TCS, 1);
TP=size(Controls_TCS, 2);
NSUB_Controls=size(Controls_TCS,3);
NSUB_AgCC=size(AgCC_TCS,3);
%create unique variable with 2 groups concatenated
All_Subjs_TCS=zeros(N_areas,TP,NSUB_Controls+NSUB_AgCC);
All_Subjs_TCS(:,:,1:NSUB_AgCC)=AgCC_TCS;%from 1 to 16
All_Subjs_TCS(:,:,NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls)=Controls_TCS;%from 17 to 44
NSUB=size(All_Subjs_TCS,3);
Ttotal=TP*NSUB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preprocessing - Remove areas with timecourses at zero

% load file containing the index of regions with timecourses at zero
% this is used to remove these regions from the subsequent analysis 
if size(Controls_TCS, 1) < 200
    load ../data/TCS/areas_zero.mat;
else
    load ../data/TCS/areas_zero_246.mat
end

All_Subjs_TCS(areas_zero,:,:) = [];

N_areas = size(All_Subjs_TCS,1);

%% 1 - Compute the Leading Eigenvectors from the BOLD datasets
disp('Processing the eigenvectors from BOLD data')
% Load here the BOLD data (which may be in different formats)
% Here the BOLD time courses in AAL parcellation are organized as cells,
% where tc_aal{1,1} corresponds to the BOLD data from subject 1 in
% condition 1 and contains a matrix with lines = N_areas and columns=Tmax.

% Parameters of the data
TR=2;  % Repetition Time (seconds)

% Preallocate variables to save FC patterns and associated information
Leading_Eig=zeros(Ttotal,N_areas); % All leading eigenvectors
Time_all=zeros(Ttotal,1); % vector with subject nr and task at each t XXXXXXXXXXX
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for s=1:NSUB
    % Get the BOLD signals from this subject in this task
    BOLD = All_Subjs_TCS(:,:,s);
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
        
        % Get the leading eigenvector
        %           [V1 ~]=eigs(iFC,1);
        
        [VV DD]=eigs(iFC,2);
        d1=DD(1,1)/sum(diag(DD));
        d2=DD(2,2)/sum(diag(DD));
        V1=d1*VV(:,1);
        V2=d2*VV(:,2);
        % Make sure the largest component is negative
        if mean(V1>0)>.5
            V1=-V1;
        elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
            V1=-V1;
        end
%         if mean(V2>0)>.5
%             V2=-V2;
%         elseif mean(V2>0)==.5 && sum(V2(V2>0))>-sum(V2(V2<0))
%             V2=-V2;
%         end
        
        % Save V1 from all frames in all fMRI sessions in Leading eig
        t_all=t_all+1; % Update time
        Time_all(t_all)=s; % belonging of each eigenvector position to subject
        % taking the first leading eigenvector
        Leading_Eig(t_all,:)=V1; %vertcat(V1,V2);
    end
end
clear BOLD tc_aal signal_filt iFC VV V1 V2 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors, using PCA

disp('Apply PCA to the eigenvectors')
% Leading_Eig is a matrix containing all the eigenvectors:
% Columns: N_areas are brain areas (variables)
% Rows: Tmax*n_Subjects are all time points (independent observations)

[U,S,V] = svd(Leading_Eig); 

%% Plot the singular values and cumulative singular values to choose the number of PCs to take
figure()
hold on
subplot(1,2,1)
semilogy(diag(S), 'k-o', 'LineWidth', 2.5)
xlabel('principal components')
title('Singular values (log scale)')
set(gca, 'FontSize', 15), axis tight, grid on

subplot(1,2,2)
plot(cumsum(diag(S))./sum(diag(S)), 'k-o', 'LineWidth', 2.5)
xlabel('principal components')
title('Cumulative singular values')
set(gca, 'FontSize', 15), axis tight, grid on 

%% recover principal components (PCs) and project data in the PCs

% determine the number of PCs that explains at least 50% of variance 
K = find(cumsum(diag(S))./sum(diag(S)) > 0.5, 1 );

Leading_Eig_PC = zeros(K, size(Leading_Eig,1));
for i=1:size(Leading_Eig,1)
    for n = 1:K
        % recover the timecourses of the PCs
        Leading_Eig_PC(n,i) = V(:,n)'*Leading_Eig(i,:)';
    end
end

%% Get the vector of groups for the K chosen PCs
V_pc = V(:,1:K);

%% Calculate probability and lifetime metrics 

P = zeros(NSUB, K); % P, Probability of Occurence
LT = zeros(NSUB, K); % LT, Mean Lifetime
for s=1:NSUB
    % every TP timepoints there is a new subject data
    % select intervals on TP timepoints 
    start_ = 1+(TP+1)*(s-1); 
    end_ = TP*s;
    [M,I] = max(Leading_Eig_PC(:,start_:end_)); 
    for k=1:K
        P(s,k) = mean(I==k);
        Time_bin=I==k;
        % Detect switches in and out of this state
        a=find(diff(Time_bin)==1);
        b=find(diff(Time_bin)==-1);
                
        % We discard the cases where state sarts or ends ON
        if length(b)>length(a)
            b(1)=[];
        elseif length(a)>length(b)
            a(end)=[];
        elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
            b(1)=[];
            a(end)=[];
        end
        if ~isempty(a) && ~isempty(b)
            durations=b-a;
        else
            durations=0;
        end
        LT(s,k)=mean(durations)*TR;
    end
end

%% find P and LT for the two groups (agcc + control)
Pctrl = P(1:NSUB_Controls,:);
Pagcc = P(NSUB_Controls+1:end,:);

LTctrl = LT(1:NSUB_Controls,:);
LTagcc = LT(NSUB_Controls+1:end,:);
%% paired t-test between two groups' probabilities and lifetimes
[sig_meas_unc_p,sig_meas_corr_p]=permutation_test([Pagcc',Pctrl']',NSUB_AgCC,0.05,19,'mean','1');
[sig_meas_unc_lt,sig_meas_corr_lt]=permutation_test([LTagcc',LTctrl']',NSUB_AgCC,0.05,19,'mean','1');

%% bar plots - probabilities of AgCC and Control
for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pagcc(:,c)); % AgCC
    Group2=squeeze(Pctrl(:,c)); % Control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'Agcc', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_meas_unc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off

    subplot(2,K,K+c)
    Group1=squeeze(LTagcc(:,c)); % AgCC
    Group2=squeeze(LTctrl(:,c)); % Control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'Agcc', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_meas_unc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])

%% Statistics between groups of AgCC: partial vs. complete 

% indeces of partial and complete AgCC
idx_partial = [1,2,10,11,12];
idx_complete = [3,4,5,6,7,8,9,13];

Ppartial=squeeze(Pagcc(idx_partial,:));
Pcomplete=squeeze(Pagcc(idx_complete,:));

LTpartial=squeeze(LTagcc(idx_partial,:));
LTcomplete=squeeze(LTagcc(idx_complete,:));

% paired t-test between the two groups' probabilities and lifetimes 
[sig_meas_unc_agcc_p,sig_meas_corr_agcc_p]=permutation_test([Ppartial',Pcomplete']',numel(idx_partial),0.05,19,'mean','1');
[sig_meas_unc_agcc_lt,sig_meas_corr_agcc_lt]=permutation_test([LTpartial',LTcomplete']',numel(idx_partial),0.05,19,'mean','1');

%% bar plots - probabilities of partial and complete AgCC
for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Ppartial(:,c)); % partial AgCC
    Group2=squeeze(Pcomplete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'pAgCC', 'cAgCC'},'FontSize', 12)
    xtickangle(-65)
    if sig_meas_unc_agcc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off

    subplot(2,K,K+c)
    Group1=squeeze(LTpartial(:,c)); % partial AgCC
    Group2=squeeze(LTcomplete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'pAgCC', 'cAgCC'},'FontSize', 12)
    xtickangle(-65)
    if sig_meas_unc_agcc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])

%% Saving
save ../results/leida/LEiDA_results_pca.mat