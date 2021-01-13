%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Plots of the FC and SC matrices to compare simulations' results 
% 
% Written by
% Ludovica Romanin - november 2020
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; 

% load Hopf model simulation results 
[file,path] = uigetfile(['../../results/hopf_model/optimizedhopf_','*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/hopf_model/',file]);

%% 
if size(FCemp,2) > 200
    vec = [1:243];
else
    vec = [1:74];
end

%% min KL, optimal G
[M,I] = min(klpstates);

%%
figure()
colormap(jet)
subplot(1,2,1)
imagesc(vec,vec,squeeze(FCsim(I,:,:)));
colorbar
title('simulated FC')

colormap(jet)
subplot(1,2,2)
imagesc(vec,vec,FCemp);
colorbar
title('empirical FC')

sgtitle('Control model - optimised for KL');

%%
% Load woCC simulation 
[file,path] = uigetfile(['../../results/hopf_model/optimizedhopf_','*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/hopf_model/',file]);

%% min KL, optimal G
[M,I] = min(klpstates);

%%
figure()
colormap(jet)
subplot(1,2,1)
imagesc(vec,vec,squeeze(FCsim(I,:,:)));
colorbar
title('simulated FC')

colormap(jet)
subplot(1,2,2)
imagesc(vec,vec,FCemp);
colorbar
title('empirical FC')

sgtitle('woCC model');

%% plot SC 

clear all; clc;
load ../../data/meanSCs/meanSC_56HC_Desikan.mat;
SC_HCP = meanSC;
load ../../data/meanSCs/meanSC_56HC_Desikan_woCC.mat;
SC_WOCC = meanSC;
%%
vec = [1:80];
figure()

colormap(jet)
subplot(1,2,1)
imagesc(vec,vec,log(SC_HCP));
colorbar
title('HCP - Control')

colormap(jet)
subplot(1,2,2)
imagesc(vec,vec,log(SC_WOCC));
colorbar
title('virtual callosotomy (woCC)')

sgtitle('Structural connectome (SC)');

%% plot BN SC matrices 

load ../../data/meanSCs/meanSC_BN_AgCC.mat;
sc_agcc = meanSC;

load ../../data/meanSCs/meanSC_BN_HC.mat;
sc_hcp = meanSC;

vec = [1:246];
figure()

colormap(jet)
subplot(1,2,1)
imagesc(vec,vec,log(sc_hcp));
colorbar
title('Control')

colormap(jet)
subplot(1,2,2)
imagesc(vec,vec,log(sc_agcc));
colorbar
title('AgCC')

sgtitle('Structural connectome (SC) - BN246');


