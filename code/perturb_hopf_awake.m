%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HOPF MODEL PERTURBATION
%
% Written by 
% Gustavo Deco
% Edited by
% Jakub Vohryzek February 2020
% Ludovica Romanin April 2020
%
% Deco, Gustavo, et al. "Awakening: Predicting external stimulation
% to force transitions between different brain states." Proceedings 
% of the National Academy of Sciences 116.36 (2019): 18088-18097.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;

addpath('functions/')

%% Loading Files
load  ../results/leida/empiricalLEiDA_K_10_DK.mat;
P1emp=mean(P1emp);
P2emp=mean(P2emp);

load  ../results/hopf_model/optimizedhopf_K_10_DK_EC.mat;

%% Perturbation Protocol
PERTURB=0:0.01:0.1;  %% This is synchronisation protocol
%PERTURB=0:-0.025:-0.4; %% This is for the noise Protocol
%% Parameters
load ../data/TCS/Controls_TCS.mat

% remove the areas with timecourses at zero 
load ../data/TCS/areas_zero.mat
Controls_TCS(areas_zero,:,:) = [];

we=0.09;
ITER=30;
N=size(Controls_TCS,1);
a=zeros(N,2);
% Coptim: from optimize_hopf_effective, it is the Effective Connectivity
% for the Global Coupling Factor (WE)
C=squeeze(Coptim(find(abs(WE-we)<0.0001),:,:));

TP=size(Controls_TCS,2);
NSUB=size(Controls_TCS,3);
TR=2;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
Tmax=NSUB*TP;
sig=0.02;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

wC = we*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

    %%%%%%%%%%%%%%% A: PERTURBATION %%%%%%%%%%%%
%%
for node=1:N/2
    node
    iwe=1;
    %% A1: PERTURBING
    for perturb=PERTURB
        a=zeros(N,2); % a will be the parameter used to add a perturbation to the Hopf model
        a(node,:)=a(node,:)+perturb;
        a(N+1-node,:)=a(N+1-node,:)+perturb;
        for iter=1:ITER
            xs=zeros(Tmax,N);
            z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            %% A2: MODEL DISCARDING: first 3000 time steps
            for t=0:dt:3000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                % Super-critical Hopf bifurcation equation
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            end
            %% A3: MODELLING
            % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
            for t=0:dt:((Tmax-1)*TR)
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                % Super-critical Hopf bifurcation equation
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
                if abs(mod(t,TR))<0.01
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            %% A4: COMPARISON
            %%%% KL dist between PTR2emp
            
            % apply LEiDA to the simulated data
            [PTRsim,Pstates]=LEiDA_fix_clusterAwakening(xs',NumClusters,Vemp,TR);
            
            % compute the KL distance between simulated and extracted states
            KLps_p(iter)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)));
        end
        KLpstatessleep_perturbed(node,iwe)=mean(KLps_p)
        iwe=iwe+1;
    end
end
%% Plotting
figure
imagesc(KLpstatessleep_perturbed);

