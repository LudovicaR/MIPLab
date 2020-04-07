%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% Written by
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
% Giulia Preti March 2020 - adaptation for AgCC analysis: 2 groups
% Ludovica Romanin March 2020 - remove small areas giving timecourses=0
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%load data matrices for the two groups
load AgCC_TCS_new.mat
load Controls_TCS.mat

%%

%parameters
N_areas=80;
TP=200;
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
 
load areas_zero.mat

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
        if mean(V2>0)>.5
            V2=-V2;
        elseif mean(V2>0)==.5 && sum(V2(V2>0))>-sum(V2(V2<0))
            V2=-V2;
        end
        
        % Save V1 from all frames in all fMRI sessions in Leading eig
        t_all=t_all+1; % Update time
        Time_all(t_all)=s; % belonging of each eigenvector position to subject
        % taking the first leading eigenvector
        Leading_Eig(t_all,:)=V1; %vertcat(V1,V2);
    end
end
clear BOLD tc_aal signal_filt iFC VV V1 V2 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors

disp('Clustering the eigenvectors into')
% Leading_Eig is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: Tmax*n_Subjects are all time points (independent observations)

% Set maximum/minimum number of clusters
% There is no fixed number of states the brain can display
% Extending depending on the hypothesis of each work

maxk=10; %%10;
mink=2; %%2;
rangeK=mink:maxk;

% Set the parameters for Kmeans clustering
Kmeans_results=cell(size(rangeK));
for k=1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' clusters'])
    [IDX, C, SUMD, D]=kmeans(Leading_Eig,rangeK(k),'Replicates',10,'MaxIter',1000,'Display','off');
    Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric column vectors
    Kmeans_results{k}.C=C;       % Cluster centroids (FC patterns)
    Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D;       % Distance from each point to every centroid
    [ss,H]=silhouette(Leading_Eig,IDX);
    sil(k)=mean(ss);
end

%% 3 - Analyse the Clustering results

% For every fMRI scan calculate probability and lifetimes of each pattern c.

% initialising P and LT

P = zeros(NSUB,maxk-mink+1,maxk);
LT = zeros(NSUB,maxk-mink+1,maxk);

for k=1:length(rangeK)

        for s=1:NSUB
            
            % Select the time points representing this subject and task
            T=(Time_all==s);
            Ctime=Kmeans_results{k}.IDX(T);
           
            for c=1:rangeK(k)
                % Calculate Probability of Occurence
                P(s,k,c)=mean(Ctime==c);
                
                % Calculate Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
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
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(s,k,c)=mean(C_Durations)*TR;
            end
        end
end

%% initialising Statistics for P and LT

P_pval=zeros(maxk-mink+1,maxk);
LT_pval=zeros(maxk-mink+1,maxk);

disp('Test significance difference between groups')

for k=1:length(rangeK)
    disp(['Now running for ' num2str(k) ' clusters'])
    for c=1:rangeK(k)
        % Compare Probabilities of Occurence
        a=squeeze(P(1:NSUB_AgCC,k,c));  % Vector containing Prob of c in AgCC
        b=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,c));  % Vector containing Prob of c in Task
        stats=permutation_htest2_np([a',b'],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        P_pval(k,c)=min(stats.pvals);
        
        % Comapre Lifetimes
        a=squeeze(LT(1:NSUB_AgCC,k,c));  % Vector containing Prob of c in Baselineline
        b=squeeze(LT(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,c));  % Vector containing Prob of c in Task
        stats=permutation_htest2_np([a',b'],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        LT_pval(k,c)=min(stats.pvals);
    end
end
disp('%%%%% LEiDA SUCCESSFULLY COMPLETED %%%%%%%')

%% 4 - Choosing FC patterns and stastistics between groups


disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])

for k=1:length(rangeK)
    Pmin_pval(k)=min(P_pval(k,1:rangeK(k)));
    %     if Pmin_pval<0.05
    %         break;
    %     end
end

[pmin_pval k]=min(Pmin_pval);
disp(['Note: The minimal K with significant difference is detected with K=' num2str(rangeK(k)) ' (p=' num2str(pmin_pval) ')'])

numsig=-1;
for k=1:length(rangeK)
    nums=length(find(P_pval(k,1:rangeK(k))<0.05));
    if nums>numsig
        numsig=nums;
    end
end
disp(['Note: The K with max number of significant difference is detected with K=' num2str(rangeK(k)) ' (n=' num2str(numsig) ')'])

disp(['Note: Silhouette optimum:'])
sil

K = input('Number of clusters: ');
Number_Clusters=K;
Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Clusters are sorted according to their probability of occurrence
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend');

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);
% Order=[1:2:N N:-2:2];

Vemp=V;

P1emp=squeeze(P(1:NSUB_AgCC,k,ind_sort));
P2emp=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort));
LT1emp=squeeze(LT(1:NSUB_AgCC,k,ind_sort));
LT2emp=squeeze(LT(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort));
PTR1emp=zeros(K);
PTR2emp=zeros(K);

n_sub1=zeros(K,1);
n_sub2=zeros(K,1);


for s=1:NSUB_AgCC
    % Select the time points representing this subject
    T=(Time_all==s);
    Ctime=Kmeans_results{k}.IDX(T);
    
    i=1;
    for c1=ind_sort
        j=1;
        for c2=ind_sort
            sumatr=0;
            for t=1:length(Ctime)-1
                if Ctime(t)==c1 && Ctime(t+1)==c2
                    sumatr=sumatr+1;
                end
            end
            if length(find(Ctime(1:length(Ctime)-1)==c1)) ~= 0
                PTR1emp(i,j)=PTR1emp(j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
            end
            j=j+1;
        end
        if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
            n_sub1=n_sub1+1;
        end
        i=i+1;
    end
end

for s=NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls
    % Select the time points representing this subject
    T=(Time_all==s);
    Ctime=Kmeans_results{k}.IDX(T);
    
    i=1;
    for c1=ind_sort
        j=1;
        for c2=ind_sort
            sumatr=0;
            for t=1:length(Ctime)-1
                if Ctime(t)==c1 && Ctime(t+1)==c2
                    sumatr=sumatr+1;
                end
            end
            if length(find(Ctime(1:length(Ctime)-1)==c1)) ~= 0
                PTR2emp(i,j)=PTR2emp(j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
            end
            j=j+1;
        end
        if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
            n_sub2=n_sub2+1;
        end
        i=i+1;
    end
end

for i=1:K
    PTR1emp(i,:)=PTR1emp(i,:)/n_sub1(i);
    PTR2emp(i,:)=PTR2emp(i,:)/n_sub2(i);
end
%% Saving

save empiricalLEiDA.mat Vemp P1emp P2emp LT1emp LT2emp PTR1emp PTR2emp Number_Clusters;

%% Plotting
disp(' ')
disp('%%% PLOTS %%%%')


figure
colormap(jet)
% Pannel A - Plot the FC patterns over the cortex
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition
% Pannel D - Plot the lifetimes of each state in each condition

for c=1:K
    subplot(5,K,c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.txt
    plot_nodes_in_cortex(V(c,1:length(labels)))
    title({['State #' num2str(c)]})
    
    subplot(5,K,K+c)
    FC_V=V(c,1:N_areas)'*V(c,1:N_areas);
    li=max(abs(FC_V(:))); % color bar simmetrica, in base al massimo o plottare nel range originale 
    imagesc(FC_V,[-li li])
    axis square
    title('FC pattern')
    ylabel('Brain area #')
    xlabel('Brain area #')
    
    subplot(5,K,2*K+c)
    Group1=squeeze(P(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'})
    if P_pval(k,ind_sort(c))<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability')
    end
    box off
    
    subplot(5,K,3*K+c)
    Group1=squeeze(P(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'})
    if LT_pval(k,ind_sort(c))<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)')
    end
    box off
end