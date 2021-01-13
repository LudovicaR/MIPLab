%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANALYSIS OF THE HOPF MODEL SIMULATION
%
% Written by
% Ludovica Romanin - june 2020
% 
% This script plots the results from the Hopf model simulation. 
% It determines the optimal G parameter according to the minimal KS and KL
% distances. 
% It is useful to assess the goodness of fitt and the metrics behavior.
% It also applies t-test to the probabilities (simulated vs. empirical) and
% plots the bar graphs. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;clc;

% load Hopf model simulation results 
[file,path] = uigetfile(['../../results/hopf_model/optimizedhopf_','*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/hopf_model/',file]);

%% discriminate between the case of a single or multiple subjects simulation
% enter 0 if only one subject was simulated 
% enter 1 for multiple subjects 

subs = str2double(input('Enter 0 if only 1 subject was simulated, 1 for multiple subjects: ','s'));
if isnan(subs) || fix(subs) ~= subs
  disp('Please enter a 0 (single subject) or 1 (multiple subjects)')
end

%% find optimal G and corresponding PMS

if subs 
    % Optimal G for KL
    [M,I] = min(klpstates_sub,[], 2);
    we_optim_kl = WE(I);
    
    % Optimal G for KS
    [M1,I1] = min(ksdist_sub, [], 2);
    we_optim_ks = WE(I1);
else
    % Optimal G for KL
    [M,I] = min(klpstates);
    we_optim_kl = WE(I);
    
    % Optimal G for KS
    [M1,I1] = min(ksdist);
    we_optim_ks = WE(I1);
end

%%
figure
hold on;
if subs 
    plot(WE,nanmean(ksdist_sub),'r');
    plot(WE,nanmean(klpstates_sub),'b');
    title('Statistics for Control model - multiple subjects')
else
    plot(WE, ksdist, 'r')
    plot(WE, klpstates, 'b')
    title('Statistics for Control model') % change title if model of woCC
end
xline(mean(we_optim_kl),'--b');
xline(mean(we_optim_ks),'--r');

xlabel('global coupling factor')
legend('KS distance', 'KL distance')

%%
figure
hold on;
if subs 
    plot(WE, nanmean(fitt_sub), 'k');
    title('Statistics for Control model - multiple subjects')
else
    plot(WE, fitt, 'k');
    title('Statistics for Control model') % change title if model of woCC
end
xline(mean(we_optim_kl),'--b');
xline(mean(we_optim_ks),'--r');

xlabel('global coupling factor')
legend('Pearson correlation', 'G_{opt}, KL', 'G_{opt}, KS', 'FontSize', 11)


%% plot when having multiple subjects
% observe distribution of optimal G parameters according to KL and KS

figure
hold on
subplot(1,2,1)
boxplot(we_optim_kl)
title('KL')
xticklabels({' '})
subplot(1,2,2)
boxplot(we_optim_ks)
title('KS')
xticklabels({' '})

%% get number of clusters
C = strsplit(file,'_');
K = str2num(C{3});
%% 
if subs
    % Optimal G, multiple subjects
    % choose KL since more stable across different subjects 
    [M,I] = min(klpstates_sub,[],2);
    we_optim_kl = WE(I);
    
    Pstatessimul = zeros(size(klpstates_sub, 1),K);
    for i=1:size(klpstates_sub, 1)
        Pstatessimul(i,:) = Pstatessimul_sub(i,I(i),:);
    end

else
    % Optimal G, one subject
    [M,I] = min(klpstates); % choose either ksdist (KS dist) or klpstates (KL dist)
    we_optim = WE(I);
end

%% load LEiDA results 

% load LEiDA results to get empirical probabilities
% make sure you load the data parcellated with the same atlas !
[file,path] = uigetfile(['../../results/leida/empiricalLEiDA_K_',num2str(K),'*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/leida/',file]);

%% t_test to assess significance of difference between PMS
if subs
    % for multiple subjects
    [sig_meas_unc_ctrl,sig_meas_corr_ctrl] = permutation_test([P2emp', Pstatessimul']',size(P2emp,1),0.05,19,'mean','2');
    [sig_meas_unc_agcc,sig_meas_corr_agcc] = permutation_test([P1emp', Pstatessimul']',size(P1emp,1),0.05,19,'mean','2');
else
    % one sample t-test
    [sig_meas_unc_ctrl,sig_meas_corr_ctrl]=permutation_onesample_ttest_2021(P2emp,Pstatessimul(I,:),0.05,19,'mean','3');
    [sig_meas_unc_agcc,sig_meas_corr_agcc]=permutation_onesample_ttest_2021(P1emp,Pstatessimul(I,:),0.05,19,'mean','3');

end

%% plot PMS, simulated and empirical

figure()
for c=1:K
    subplot(1,K,c) % 'position', [left bottom width height]
    if subs
        Group1=squeeze(Pstatessimul(:,c)); % simulated
    else
        Group1=squeeze(Pstatessimul(I,c)); % simulated
    end
    Group2=squeeze(P2emp(:,c)); % empirical control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    if c==8
        title({['Simulated Control vs. Empirical Control'] ['State #' num2str(c)]}, 'FontSize', 12);
    else
        title({['State #' num2str(c)]}, 'FontSize', 12);
    end
    
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'sim', 'emp'},'FontSize', 12)
    xtickangle(-60)
    if sig_meas_unc_ctrl(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_meas_corr_ctrl(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    %box off
end

figure
for c=1:K
    subplot(1,K,c)
    if subs
        Group1=squeeze(Pstatessimul(:,c)); % simulated
    else
        Group1=squeeze(Pstatessimul(I,c)); % simulated 
    end
    Group2=squeeze(P1emp(:,c)); % empirical control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    if c==8
        title({['Simulated Control vs. Empirical AgCC'] ['State #' num2str(c)]}, 'FontSize', 12);
    else
        title({['State #' num2str(c)]}, 'FontSize', 12);
    end
    
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'sim', 'emp'},'FontSize', 12)
    xtickangle(-60)
    if sig_meas_unc_agcc(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_meas_corr_agcc(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
end

