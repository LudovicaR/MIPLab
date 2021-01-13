%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HOPF MODEL OPTIMISATION - grid search 
%
% Written by 
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
% Ludovica Romanin March 2020
% Ludovica Romanin November 2020: grid search for G and a
%
% Deco, Gustavo, et al. "Awakening: Predicting external stimulation
% to force transitions between different brain states." Proceedings 
% of the National Academy of Sciences 116.36 (2019): 18088-18097.
%
% Hopf model simulation (see optimize_hopf_noEC.m for more details)
% The result is used to explore the fitting measures, like in the paper:
% (Deco, 2017) The dynamics of resting fluctuations
%    (doi:10.1038/s41598-017-03073-5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
clear all;clc;
tic

addpath('functions/')

%% Load LEiDA results

% load file with LEiDA results
% load empiricalLEiDA_*.mat file
load  ../results/leida/empiricalLEiDA_.mat; 

P2emp=mean(P2emp); % mean probability of occurrence of the Control group

%% Load empirical data 
load ../data/TCS/Controls_TCS.mat;

%% Load Structural Connectivity
% The structural connectivity matrix is obtained using diffusion MRI and
% tractography. 
% The number of coupled dynamical units is equal to the number of cortical 
% and subcortical areas from the AAL (atlas) parcellation.

% choose SC matrix to test for the model simulation
load ../data/meanSCs/meanSC_56HC_Desikan.mat; 
C=meanSC;
C=C/max(max(C))*0.2; 

%% remove the areas with timecourses at zero from the SC matrix

% load file containing the index of regions with timecourses at zero
% this is used to remove these regions from the subsequent analysis 
if size(Controls_TCS, 1) < 200
    load ../data/TCS/areas_zero.mat;
else
    load ../data/TCS/areas_zero_246.mat
end

C(areas_zero,:) = [];
C(:,areas_zero) = [];

% remove the areas with timecourses at zero 
Controls_TCS(areas_zero,:,:) = [];

%% Initialize parameters
% to reduce the computational cost, take only 5 subjects from empirical data
irand = randi([1 28],1,5); % to get 10 random subject data
Controls_TCS = Controls_TCS(:,:,irand);

TP=size(Controls_TCS,2);
NSUB=size(Controls_TCS,3);
TR=2;  % Repetition Time (seconds)
NumClusters=Number_Clusters; % parameter stored in empiricalLEiDA.mat

%% Filtering
delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;            % lowpass frequency of filter
fhi = fnq-0.001;      % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);   % construct the filter

flp = .04;            % lowpass frequency of filter
fhi = .07;            % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

    %%%%%%%%%%%%%%% A: EXPERIMENTAL %%%%%%%%%%%%
%%
% Extracting FC FCD and metastability of data
% FCD: Functional Connectivity Dynamics 
kk=1;
insub=1;

N=size(Controls_TCS,1);
% tril: extract lower triangular part; 
Isubdiag = find(tril(ones(N),-1));

for nsub=1:NSUB
    signaldata=Controls_TCS(:,:,nsub);
    Phase_BOLD_data=zeros(N,TP);
    timeseriedata = zeros(N,TP);
     %% A1a: ORDER PARAMETER and iFC

    for seed=1:N
        x=demean(detrend(signaldata(seed,:)));
        x(find(x>3*std(x)))=3*std(x);
        x(find(x<-3*std(x)))=-3*std(x);
        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
        Phase_BOLD_data(seed,:) = angle(hilbert(timeseriedata(seed,:)));
    end
    T=10:TP-10;
    for t=T
        kudata=sum(complex(cos(Phase_BOLD_data(:,t)),sin(Phase_BOLD_data(:,t))))/N;
        syncdata(t-9)=abs(kudata);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD_data(i,t),Phase_BOLD_data(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastabilitydata2(nsub)=std(syncdata);
    %% A1b: DYNAMIC FUNCTIONAL CONNECTIVITY (DFC)
    
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            % phfcddata: FCD for empirical data
            phfcddata(kk)=dot(p1,p2)/norm(p1)/norm(p2); % cosine similary for KS distance
            kk=kk+1;
        end
    end
    
    for t=1:TP
        for n=1:N
            for p=1:N
                % iFC: Instantaneous Functional Connectivity
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end
    FCphasesemp2(nsub,:,:)=squeeze(mean(iFC));
end
FCphasesemp=squeeze(mean(FCphasesemp2));
metastabilitydata=mean(metastabilitydata2);

%% A2: SPECTRAL ANALYSIS
% Extracting peak of data power spectra for determining omega (Hopf)
for nsub=1:NSUB
    clear PowSpect PowSpect2;
    N = size(Controls_TCS, 1);
    T = size(Controls_TCS, 2);
    Isubdiag = find(tril(ones(N),-1));
    Ts = T*TR;
    freq = (0:T/2-1)/Ts;
    signaldata=Controls_TCS(:,:,nsub);
    FCemp2(nsub,:,:)=corrcoef(signaldata');
    
    %%%%
    
    [aux minfreq]=min(abs(freq-0.04)); % min frequency of the cut-off
    [aux maxfreq]=min(abs(freq-0.07)); % max frequency of the cut-off
    nfreqs=length(freq);
    
    for seed=1:N
        % demean: function of the project, removes the mean value 
        % detrend: matlab tool, removes the best straight-line fit linear trend from the data
        x=detrend(demean(signaldata(seed,:))); % demeaning and detrending
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        % fft: compute the Discrete Fourier Transform 
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(T/2)).^2/(T/TR);
        ts2 =zscore(filtfilt(bfilt,afilt,x));
        pw2 = abs(fft(ts2));
        PowSpect2(:,seed,insub) = pw2(1:floor(T/2)).^2/(T/TR);
    end
    insub=insub+1;
end

Power_Areas=mean(PowSpect,3);
Power_Areas2=mean(PowSpect2,3);
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
    Power_Areas2(:,seed)=gaussfilt(freq,Power_Areas2(:,seed)',0.01);
    vsig(seed)=sum(Power_Areas2(minfreq:maxfreq,seed))/sum(Power_Areas2(:,seed));
end

vmax=max(vsig);
vmin=min(vsig);

[maxpowdata,index]=max(Power_Areas); % maximum power for a regions at a given frequency
f_diff = freq(index); % the frequency (the intrinsic frequency for each region)
FCemp=squeeze(mean(FCemp2));

clear PowSpect PowSpect2 Power_Areas Power_Areas2;

%% Modelling

omega = repmat(2*pi*f_diff',1,2); % multiplying by 2pi and setting two column for x and y modelling eqn
omega(:,1) = -omega(:,1); % setting the first column negative as omega is negative in x-eqn

dt=0.1*TR/2;
Tmax=TP*NSUB;
sig=0.02; % standard deviation of the noise term (corresponds to beta in the equation)
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

                  %%%%%%%%%%%% B: SIMULATION %%%%%%%%%%%%
%%
ia=1; % count for the bifurcation param

WE=0:0.01:1; % 0:0.005:0.15;  %% Global Coupling Factor G

a_coeff = [-0.5:0.02:0.5]; % test different bifurcation parameters (as in Jobst et al., 2017)
a=ones(N,2); % bifurcation parameter vector 

Na = length(a_coeff);
NWE=length(WE);
PTRsimul=zeros(Na, NWE,NumClusters,NumClusters);
Pstatessimul=zeros(Na, NWE,NumClusters);

FCsim = zeros(Na, NWE, N, N);

for a_i=a_coeff
    iwe=1; % count for the global wieght iteration
    for we=WE % loops over changing coupling constant G
    %% B2: FINAL SIMULATION
    
    wC = we*C; % multiplies SC with we(constant G)
    sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
    
    xs=zeros(Tmax,N);
    % number of iterations, 100
    z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    %% B2a: MODEL DISCARD: first 3000 time steps
    for t=0:dt:3000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a_i*a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    %% B2b: MODEL 
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a_i*a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    %% B2c: COMPUTE COMPARISON BETWEEN STATIC EMPIRICAL AND SIMULATED  FC

    FC_simul=corrcoef(xs(1:nn,:));
    FCsim(ia, iwe, :, :) = FC_simul; 
    cc = corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)));
    fitt(ia, iwe)=cc(2);
    %% B2d: COMPUTE: HILBERT TRANSFORM

    BOLD=xs';
    Phase_BOLD=zeros(N,nn);
    signal_filt=zeros(N,nn);
    for seed=1:N
        BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
    end
    %% B2e: ORDER PARAMETER

    T=10:Tmax-10;
    for t=T
        % ku = sum of all elements in complex array, the cos and sin of
        % BOLD phases, divided by the number of signals; 
        % similar to an average
        ku=sum(complex(cos(Phase_BOLD(:,t)),sin(Phase_BOLD(:,t))))/N;
        sync(t-9)=abs(ku);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD(i,t),Phase_BOLD(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    % Metastability: the quality of systems to temporarily persist in an 
    % existing equilibrium despite slight perturbations.
    metastability(ia, iwe)=abs(metastabilitydata-std(sync));
    %% B2f: DYNAMIC FUNCTIONAL CONNECTIVITY (DFC)

    kk=1;
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            % phfcd is the FCD for simulated data 
            phfcd(kk)=dot(p1,p2)/norm(p1)/norm(p2); % cosine  similarity for KS distance
            kk=kk+1;
        end
    end
    
    % kstest2: two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis
    % test. It determines if independent random samples (inputs) are drawn 
    % from the same underlying continuous population. 
    % H is the result of the hypothesis test (=0 or =1 for rejection).
    % 'ksdist' is the KS test statistic defined according to the TAIL.
    % The test statistic is equal to max|S1(x) - S2(x)| for TAIL='unequal', 
    % max[S1(x) - S2(x)] for TAIL='larger', and max[S2(x) - S1(x)] for 
    % TAIL='smaller'.
    % P is the P-value, compared to the significance level of 5%.
    % pdfcd is for the model (simulated)
    % phfcddata is for the empirical 
    [H,P,ksdist(ia, iwe)]=kstest2(phfcd,phfcddata);
    
                          %%%%%%%%%%%%  C: COMPARISON %%%%%%%%%%%%
    %% C1a: PROBABILISTIC STATE SPACE

    [PTRsim,Pstates]=LEiDA_fix_clusterAwakening(xs',NumClusters,Vemp,TR);  % Vemp, parameter from empiricalLEiDA.mat
    
    %% C1b: KL-DISTANCE BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF OCCURENCE
    klpstates(ia, iwe)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
    
    %% C1c:EXTRA FITTING BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF TRANSITION
    
    % to test the significance of the two metrics for the model fitting 
    
    kldist_pms(ia, iwe)=KLdist(PTR2emp,PTRsim); % PTR2emp, parameter from empiricalLEiDA.mat
    entropydist_pms(ia, iwe)=EntropyMarkov(PTR2emp,PTRsim); 
    
    PTRsimul(ia, iwe,:,:)=PTRsim;
    %%%
    Pstatessimul(ia, iwe,:)=Pstates;
    
    iwe=iwe+1;
    
    ksdist
    klpstates

    end
    ia=ia+1;
end


%% Saving

save('../results/gridSearch/HopfModel_results_grid_.mat', '-v7.3');

elapsed_time = toc;
toc
