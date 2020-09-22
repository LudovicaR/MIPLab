clear all; 
load LEiDA_results.mat;

%% For plots, put areas removed back to zero
V = Vemp;

% set areas that were removed back to zero for the plots
zero_vec = zeros(size(All_Subjs_TCS(1,:,:)));
for i=1:length(areas_zero)
    area = areas_zero(i);
    All_Subjs_TCS = cat(1,All_Subjs_TCS(1:area-1,:,:),zero_vec,All_Subjs_TCS(area:end,:,:));
    V = cat(2,V(:,1:area-1),zeros(1,size(V,1))',V(:,area:end));
end

N_areas = size(All_Subjs_TCS,1);

% put labels back at 80
labels = 1:80;

%% plot of LEiDA results: FC connections, FC matrices, Probabilites and Lifetimes
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
    plot_nodes_in_cortex(V(c,1:length(labels)),0)
    title({['State #' num2str(c)]}, 'FontSize', 12)
    
    subplot(5,K,K+c)
    FC_V=V(c,1:N_areas)'*V(c,1:N_areas);
    li=max(abs(FC_V(:))); % color bar simmetrica, in base al massimo o plottare nel range originale 
    imagesc(FC_V,[-li li])
    axis square
    title('FC pattern','FontSize', 12)
    ylabel('Brain area #','FontSize', 12)
    xlabel('Brain area #','FontSize', 12)
    
    subplot(5,K,2*K+c)
    Group1=squeeze(P(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'},'FontSize', 12)
    if P_pval(k,ind_sort(c))<0.05/K
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
    
    subplot(5,K,3*K+c)
    Group1=squeeze(LT(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(LT(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'},'FontSize', 12)
    if LT_pval(k,ind_sort(c))<0.05/K
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

%% plots of relevant states FC connections (view of three planes)

figure
colormap(jet)
states_list = [3, 9, 10];
for c=1:3
    subplot(3,3,c) 
    plot_nodes_in_cortex(V(states_list(c),1:length(labels)),0)
    title({['State #' num2str(states_list(c))] [ ] ['axial']},'FontSize', 12) % dorsal
end

for c=4:6
    subplot(3,3,c)
    plot_nodes_in_cortex(V(states_list(c-3),1:length(labels)),1)
    title({['sagittal (L)']},'FontSize', 12)  % lateral (L)
end

for c=7:9
    subplot(3,3,c)
    plot_nodes_in_cortex(V(states_list(c-6),1:length(labels)),2)
    title({['coronal']},'FontSize', 12) % front
end

%% plot significant state to detect regions involved 

state = 9;

figure
colormap(jet)

V_state = V(state,1:length(labels));

subplot(3,3,1) 
plot_nodes_in_cortex(V(state,1:length(labels)),0)
title({['axial']},'FontSize', 12) % dorsal

subplot(3,3,2)
plot_nodes_in_cortex(V(state,1:length(labels)),1)
title({['State #' num2str(state)] [ ] ['sagittal (L)']},'FontSize', 12)  % lateral (L)

subplot(3,3,3)
plot_nodes_in_cortex(V(state,1:length(labels)),2)
title({['coronal']},'FontSize', 12) % front

%% number of labels of the regions involved in state 9
labels(find(V_state>0))

%% test significance of complete vs. partial AgCC

% indeces of partial and complete AgCC
idx_partial = [1,2,10,11,12];
idx_complete = [3,4,5,6,7,8,9,13];

Pval_agcc_grps=zeros(1,k);

for c=1:rangeK(k)
        % Compare Probabilities of Occurence
        prt=squeeze(Pemp_partial(:,c));  % Vector containing Prob of c in AgCC_complete
        cmp=squeeze(Pemp_complete(:,c));  % Vector containing Prob of c in AgCC_partial
        stats=permutation_htest2_np([prt',cmp'],[ones(1,numel(prt)) 2*ones(1,numel(cmp))],1000,0.05,'ttest');
        Pval_agcc_grps(c)=min(stats.pvals);
end

for c=1:K
    subplot(1,K,c)
    Group1=squeeze(Pemp_partial(:,c)); % partial AgCC
    Group2=squeeze(Pemp_complete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'prt', 'cmp'},'FontSize', 12)
    if Pval_agcc_grps(c)<0.05/K
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
end
