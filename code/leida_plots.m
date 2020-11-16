clear all; 
load("../results/LEiDA_results_K_16_N800.mat");

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
labels = 1:819;

%% regions involved in the different states

atlas = importdata('dbs80symm_labels.txt'); 
regions_pos = {};
regions_neg = {};

for i=1:K
    regions_pos{i} = atlas(labels(find(V(i,:)>0)));
    regions_neg{i} = atlas(labels(find(V(i,:)<0)));
end

%% regions at zero
regions_zero = atlas(labels(find(V(1,:)==0)));

%% test significance of complete vs. partial AgCC

% indeces of partial and complete AgCC
idx_partial = [1,2,10,11,12];
idx_complete = [3,4,5,6,7,8,9,13];

LTemp_partial=squeeze(LT1emp(idx_partial,:));
LTemp_complete=squeeze(LT1emp(idx_complete,:));

[sig_meas_unc_agcc_p,sig_meas_corr_agcc_p]=permutation_test([Pemp_partial',Pemp_complete']',numel(idx_partial),0.05,999,'mean','1');
[sig_meas_unc_agcc_lt,sig_meas_corr_agcc_lt]=permutation_test([LTemp_partial',LTemp_complete']',numel(idx_partial),0.05,999,'mean','1');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_partial(:,c)); % partial AgCC
    Group2=squeeze(Pemp_complete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'prt', 'cmp'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_agcc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off

    subplot(2,K,K+c)
    Group1=squeeze(LTemp_partial(:,c)); % partial AgCC
    Group2=squeeze(LTemp_complete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'prt', 'cmp'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_agcc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

%% test significance of AgCC vs. Control 

[sig_meas_unc_p,sig_meas_corr_p]=permutation_test([P1emp', P2emp']',numel(P1emp(:,1)),0.05,999,'mean','1');
[sig_meas_unc_lt,sig_meas_corr_lt]=permutation_test([LT1emp', LT2emp']',numel(LT1emp(:,1)),0.05,999,'mean','1');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(P(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
    
    subplot(2,K,c+K)
    Group1=squeeze(LT(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(LT(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

%% test significance of AgCC complete vs. Control
[sig_meas_unc_p,sig_meas_corr_p]=permutation_test([Pemp_complete', P2emp']',numel(Pemp_complete(:,1)),0.05,999,'mean','1');
[sig_meas_unc_lt,sig_meas_corr_lt]=permutation_test([LTemp_complete', LT2emp']',numel(LTemp_complete(:,1)),0.05,999,'mean','1');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_complete(:,c)); % AgCC complete
    Group2=squeeze(P2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC cmp', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
    
    subplot(2,K,K+c)
    Group1=squeeze(LTemp_complete(:,c)); % AgCC complete
    Group2=squeeze(LT2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC cmp', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

%% test significance of AgCC partial vs. Control
[sig_meas_unc_p,sig_meas_corr_p]=permutation_test([Pemp_partial', P2emp']',numel(Pemp_partial(:,1)),0.05,999,'mean','1');
[sig_meas_unc_lt,sig_meas_corr_lt]=permutation_test([LTemp_partial', LT2emp']',numel(LTemp_partial(:,1)),0.05,999,'mean','1');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_partial(:,c)); % AgCC complete
    Group2=squeeze(P2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC prt', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
    
    subplot(2,K,K+c)
    Group1=squeeze(LTemp_partial(:,c)); % AgCC complete
    Group2=squeeze(LT2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC prt', 'Controls'},'FontSize', 12)
    xtickangle(-50)
    if sig_meas_unc_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

%% functional connectivity matrices for each states
for c=1:K
    if K < 9
        subplot(2,8,c)
        FC_V=V(c,1:N_areas)'*V(c,1:N_areas);
        imagesc(FC_V)
        axis square
        title({['State #' num2str(c)]}, 'FontSize', 12)
        ylabel('Brain area #')
        xlabel('Brain area #')
    else
        subplot(2,8,c)
        FC_V=V(c,1:N_areas)'*V(c,1:N_areas);
        imagesc(FC_V)
        axis square
        title({['State #' num2str(c)]}, 'FontSize', 12)
        ylabel('Brain area #')
        xlabel('Brain area #')
    end
    
end

