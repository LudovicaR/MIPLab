%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HOPF MODEL OPTIMISATION
%
% Written by 
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
%
% Deco, Gustavo, et al. "Awakening: Predicting external stimulation
% to force transitions between different brain states." Proceedings 
% of the National Academy of Sciences 116.36 (2019): 18088-18097.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
clear all;

%% Load LEiDA results
load  empiricalLEiDA.mat;

P1emp=mean(P1emp);
P2emp=mean(P2emp);

%% Load Structural Connectivity
load sc90.mat;
C=sc90;
C=C/max(max(C))*0.2;

%% Load data
load data_Awake.mat;

TSmax=240;
NSUB=18;
TR=2.08;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

%% Filtering
delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;           % lowpass frequency of filter
fhi = fnq-0.001;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);   % construct the filter

flp = .04;           % lowpass frequency of filter
fhi = .07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for nsub=1:NSUB
    [N, Tmax0]=size(X{1,nsub});
    if Tmax0<TSmax
        X{1,nsub}=[];
    end
end
X=X(~cellfun('isempty',X));
n_Subjects=size(X,2);


    %%%%%%%%%%%%%%% A: EXPERIMENTAL %%%%%%%%%%%%
%%
% Extracting FC FCD and metastability of data
kk=1;
insub=1;
TSmax=1000;
N=90;
Isubdiag = find(tril(ones(N),-1));
Tmaxtotal=0;
for nsub=1:n_Subjects
    [N, Tmax0]=size(X{1,nsub});
    Tmax=min(TSmax,Tmax0);
    Tmaxtotal=Tmaxtotal+Tmax;
    signaldata = X{1,nsub};
    signaldata=signaldata(:,1:Tmax);
    Phase_BOLD_data=zeros(N,Tmax);
    timeseriedata=zeros(N,Tmax);
     %% A1a: ORDER PARAMETER and iFC

    for seed=1:N
        x=demean(detrend(signaldata(seed,:)));
        x(find(x>3*std(x)))=3*std(x);
        x(find(x<-3*std(x)))=-3*std(x);
        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
        Phase_BOLD_data(seed,:) = angle(hilbert(timeseriedata(seed,:)));
    end
    T=10:Tmax-10;
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
            phfcddata(kk)=dot(p1,p2)/norm(p1)/norm(p2);
            kk=kk+1;
        end
    end
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
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
for nsub=1:n_Subjects
    clear PowSpect PowSpect2;
    [N, Tmax0]=size(X{1,nsub});
    Isubdiag = find(tril(ones(N),-1));
    Tmax=min(TSmax,Tmax0);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    signaldata = X{1,nsub};
    signaldata=signaldata(:,1:Tmax);
    FCemp2(nsub,:,:)=corrcoef(signaldata');
    
    %%%%
    
    [aux minfreq]=min(abs(freq-0.04)); % min frequency of the cut-ff
    [aux maxfreq]=min(abs(freq-0.07)); % max frequency of the cut-ff
    nfreqs=length(freq);
    
    
    for seed=1:N
        x=detrend(demean(signaldata(seed,:))); % demeaning and detrending
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(TT/2)).^2/(TT/TR);
        ts2 =zscore(filtfilt(bfilt,afilt,x));
        pw2 = abs(fft(ts2));
        PowSpect2(:,seed,insub) = pw2(1:floor(TT/2)).^2/(TT/TR);
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
Tmax=TSmax*n_Subjects;
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

for we=WE % loops over changing coupling constant G
    minm=100;
    Cnew=C;
    %% B1: Optimize - ANEC
    for iter=1:150 % iterations for the ANEC optimisation
        wC = we*Cnew; % multiplies SC with we(constant G)
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        xs=zeros(Tmax,N);
        % number of iterations, 100
        z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        %% B1a: MODEL DISCARD: first 3000 time steps
        for t=0:dt:3000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        end
        %% B1b: MODEL:
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
        %% B1c: COMPUTE: hilbert transform
        BOLD=xs';
        signal_filt=zeros(N,nn);
        Phase_BOLD=zeros(N,nn);
        for seed=1:N
            BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
        end
        %% B1d: COMPUTE: iFC
        for t=1:nn
            for n=1:N
                for p=1:N
                    iFC(t,n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
        end
        FCphases=squeeze(mean(iFC));
        
        %% B1e: COMPUTE: ANEC
        % Update Effective Connectivity Matrix Cnew
        for i=1:N
            for j=i+1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+0.01*(FCphasesemp(i,j)-FCphases(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                    Cnew(j,i)=Cnew(i,j);
                end
            end
        end
        
        Cnew=Cnew/max(max(Cnew))*0.2;
        
        D = abs(FCphasesemp-FCphases).^2;
        MSE = sum(D(:))/numel(FCphases);
        if MSE<0.01
            break;
        end
        
    end 
    Coptim(iwe,:,:)=Cnew;  %% Effective Connectivity for G (we)

    %% B2: FINAL SIMULATION
    
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
        ku=sum(complex(cos(Phase_BOLD(:,t)),sin(Phase_BOLD(:,t))))/N;
        sync(t-9)=abs(ku);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD(i,t),Phase_BOLD(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastability(iwe)=abs(metastabilitydata-std(sync));
    %% B2f: DYNAMIC FUNCTIONAL CONNECTIVITY (DFC)

    kk=1;
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            phfcd(kk)=dot(p1,p2)/norm(p1)/norm(p2);
            kk=kk+1;
        end
    end
    
    [H,P,ksdist(iwe)]=kstest2(phfcd,phfcddata);
    
                          %%%%%%%%%%%%  C: COMPARISON %%%%%%%%%%%%
    %% C1a: PROBABILISTIC STATE SPACE

    [PTRsim,Pstates]=LEiDA_fix_clusterAwakening(xs',NumClusters,Vemp,TR);   
    
    %% C1b: KL-DISTANCE BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF OCCURENCE
    klpstatessleep(iwe)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
    klpstatesawake(iwe)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)));
    
    %% C1c:EXTRA FITTING BETWEEN EMPIRICAL AND SIMULATED PROBABILITY OF TRANSITION
    
    kldistsleep(iwe)=KLdist(PTR2emp,PTRsim);
    kldistawake(iwe)=KLdist(PTR1emp,PTRsim);
    entropydistsleep(iwe)=EntropyMarkov(PTR2emp,PTRsim);
    entropydistawake(iwe)=EntropyMarkov(PTR1emp,PTRsim);
    
    PTRsimul(iwe,:,:)=PTRsim;
    %%%
    Pstatessimul(iwe,:)=Pstates;
    
    iwe=iwe+1;
    
    ksdist
    klpstatesawake

end
%% Saving
save optimizedhopfawake.mat WE PTRsimul Pstatessimul metastability ksdist klpstatessleep klpstatesawake kldistsleep kldistawake entropydistawake entropydistsleep fitt Coptim n_Subjects f_diff;

                    %%%%%%%%%%%%  D: VISUALISATION %%%%%%%%%%%%
%% D1: PLOTTING

figure
plot(WE,fitt,'b');
hold on;
plot(WE,kldistawake,'k'); %% extra
plot(WE,entropydistawake,'k'); %%   extra
figure
plot(WE,metastability,'r');
hold on;
plot(WE,ksdist,'c');
plot(WE,klpstatesawake,'k');
plot(WE,klpstatessleep,'b');

% 
% figure
% plot(WE,fitt,'b');
% figure
% plot(WE,kldistsleep,'r');
% hold on;
% plot(WE,kldistawake,'k');
% figure
% plot(WE,entropydistsleep,'r');
% hold on;
% plot(WE,entropydistawake,'k');    
% figure
% plot(WE,klpstatessleep,'r');
% hold on;
% plot(WE,klpstatesawake,'k');   
        