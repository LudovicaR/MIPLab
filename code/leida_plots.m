
load LEiDA_results.mat;

%% For plots, put areas removed back to zero
V = Vemp;

% set areas that were removed back to zero for the plots
zero_vec = zeros(size(All_Subjs_TCS(1,:,:)));
for i=1:length(areas_zero)
    area = areas_zero(i);
    All_Subjs_TCS = cat(1,All_Subjs_TCS(1:area,:,:),zero_vec,All_Subjs_TCS(area+1:end,:,:));
    V = cat(2,V(:,1:area),zeros(1,size(V,1))',V(:,area+1:end));
end

N_areas = size(All_Subjs_TCS,1);

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
    Group1=squeeze(LT(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(LT(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
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

%% plots of relevant states FC connections (view of three planes)

figure
colormap(jet)
states_list = [3, 9, 10];
for c=1:3
    subplot(3,3,c) 
    plot_nodes_in_cortex(V(states_list(c),1:length(labels)),0)
    title({['State #' num2str(states_list(c))] [ ] ['dorsal']},'FontSize', 12) % axial plane
end

for c=4:6
    subplot(3,3,c)
    plot_nodes_in_cortex(V(states_list(c-3),1:length(labels)),1)
    title({['lateral (L)']},'FontSize', 12)  % sagittal plane
end

for c=7:9
    subplot(3,3,c)
    plot_nodes_in_cortex(V(states_list(c-6),1:length(labels)),2)
    title({['front']},'FontSize', 12) % coronal plane
end

%% compare complete vs. partial AgCC

% indeces of partial and complete AgCC
idx_partial = [1,2,10,11,12];
idx_complete = [3,4,5,6,7,8,9,13];

n_partial = size(idx_partial,2);
n_complete = size(idx_complete,2);

figure

for c = 1:K
    
    subplot(2,K,c)
    Group1=squeeze(P(idx_partial,k,ind_sort(c))); % AgCC partial
    Group2=squeeze(P(idx_complete,k,ind_sort(c))); % AgCC complete
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'partial', 'complete'})
    if P_pval(k,ind_sort(c))<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('AgCC probability')
    end
    box off
    title({['State #' num2str(c)]})
    
    subplot(2,K,K+c)
    Group1=squeeze(LT(idx_partial,k,ind_sort(c))); % AgCC partial
    Group2=squeeze(LT(idx_complete,k,ind_sort(c))); % AgCC complete
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'partial', 'complete'})
    if LT_pval(k,ind_sort(c))<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('AgCC Lifetime (seconds)')
    end
    box off
end



