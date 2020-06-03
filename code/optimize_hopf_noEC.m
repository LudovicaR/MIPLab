%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HOPF MODEL OPTIMISATION
%
% Written by 
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
% Ludovica Romanin March 2020
%
% Deco, Gustavo, et al. "Awakening: Predicting external stimulation
% to force transitions between different brain states." Proceedings 
% of the National Academy of Sciences 116.36 (2019): 18088-18097.
%
% OUTPUTs: 
% - optimizedhopfawake.mat: mat file with all the outputs listed below.
% - WE: global coupling factor (G)
% - PTRsimul: 
%       output of LEiDA_fix_clusterAwakening.m, the Probability of transition
% - Pstatessimul: 
%       output of LEiDA_fix_clusterAwakening.m, the Probability states
% - metastability: 
%       the standard deviation of the sum of complex array composed of cos 
%       and sin of the phase of BOLD signals.
% - ksdist: the Kolmogorov-Smirnov distance statistic.
% - klpstates: symmetrized K-L distance between empirical and simulated in control group.
% - kldist_pms: symmetrized K-L for the transition state  
% - entropydist_pms: Markov Entropy 
% - fitt: 
%       correlation coefficient between empirical and simulated Functional Connectivites (FC)
% - n_Subjects: the number of subjects. 
% - f_diff: the intrinsic frequency for each region
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
clear all;

%% Load LEiDA results
load  empiricalLEiDA.mat; % file with results computed in LEiDA_2Conditions.m

% parameters from empiricalLEiDA.mat
P2emp=mean(P2emp); % mean probability of occurrence of the Control group

%% Load Structural Connectivity
% The structural connectivity matrix is obtained using diffusion MRI and
% tractography. 
% The number of coupled dynamical units is equal to the number of cortical 
% and subcortical areas from the AAL (atlas) parcellation.
load meanSC_56HC_Desikan_woCC.mat;
C=meanSC;
C=C/max(max(C))*0.2;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
C(areas_zero,:) = [];
C(:,areas_zero) = [];

%% Load data

load Controls_TCS.mat;
% remove the areas with timecourses at zero 
Controls_TCS(areas_zero,:,:) = [];

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

phfcddata_emp = zeros(1, (TP-21)*NSUB*((TP-21-1)/2));

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

phfcddata_emp = phfcddata;

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
sig=0.02; % standard deviation of the noise term
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

                  %%%%%%%%%%%% B: SIMULATION %%%%%%%%%%%%
%%
iwe=1; % count for the global wieght iteration
WE=0:0.01:1; % 0:0.005:0.15;  %% Global Coupling Factor G
a=zeros(N,2); % a bifurcation parameter if = 0 => at the bifurcation

NWE=length(WE);
PTRsimul=zeros(NWE,NumClusters,NumClusters);
Pstatessimul=zeros(NWE,NumClusters);

pat_size = N*(N-1) + (TP-21) - 2;

phfcd_sim = zeros(NWE, 1, pat_size*(pat_size-1)/2);

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
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    %% B2b: MODEL 
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    %% B2c: COMPUTE COMPARISON BETWEEN STATIC EMPIRICAL AND SIMULATED  FC

    FC_simul=corrcoef(xs(1:nn,:));
    cc = corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)));
    fitt(iwe)=cc(2);
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
    metastability(iwe)=abs(metastabilitydata-std(sync));
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
    [H,P,ksdist(iwe)]=kstest2(phfcd,phfcddata);
    
    phfcd_sim(iwe,:,:) = phfcd;
    
                          %%%%%%%%%%%%  C: COMPARISON %%%%%%%%%%%%
    %% C1a: PROBABILISTIC STATE SPACE

    [PTRsim,Pstates]=LEiDA_fix_clusterAwakening(xs',NumClusters,Vemp,TR);  % Vemp, parameter from empiricalLEiDA.mat
    
    %% C1b: KL-DISTANCE BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF OCCURENCE
    klpstates(iwe)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
    
    %% C1c:EXTRA FITTING BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF TRANSITION
    
    % to test the significance of the two metrics for the model fitting 
    
    kldist_pms(iwe)=KLdist(PTR2emp,PTRsim); % PTR2emp, parameter from empiricalLEiDA.mat
    entropydist_pms(iwe)=EntropyMarkov(PTR2emp,PTRsim); 
    
    PTRsimul(iwe,:,:)=PTRsim;
    %%%
    Pstatessimul(iwe,:)=Pstates;
    
    iwe=iwe+1;
    
    ksdist
    klpstates

end
%% Saving
save optimizedhopfawake_56HC_woCC.mat WE PTRsimul Pstatessimul metastability ksdist klpstates kldist_pms entropydist_pms fitt NSUB f_diff phfcddata_emp;

save('fcd_result_56HC_woCC.mat', 'phfcd_sim', '-v7.3');

save HopfModel_results_56HC_woCC.mat

                    %%%%%%%%%%%%  D: VISUALISATION %%%%%%%%%%%%
%% D1: PLOTTING

% WE: global coupling factor
% fitt: correlation coefficient between empirical and simulated Functional Connectivites (FC)
% kldistawake: symmetrized K-L for the transition state
% entropydistawake: Markov Entropy for the state. 
% metastability: the quality of systems to temporarily persist in an existing equilibrium despite slight perturbations.
% ksdist: Kolmogorov-Smirnov distance statistic
% klpstatesawake: symmetrized K-L distance between empirical and simulated
