clear all;
load meanSC_56HC_Desikan_woCC.mat;
C=meanSC;

% remove the areas with timecourses at zero from the SC matrix
load areas_zero.mat
C(areas_zero,:) = [];
C(:,areas_zero) = [];

load  optimizedhopfawake_56HC_woCC.mat;
load empiricalLEiDA.mat P1emp P2emp PTR1emp PTR2emp;

% Optimal G for KL
[M,I] = min(klpstates);
we_optim_kl = WE(I);

% Optimal G for KS
[M1,I1] = min(ksdist);
we_optim_ks = WE(I1);

%% t_test to assess significance of difference between PMS

K = 10;
Pval_ctrl = zeros(1,K);
Pval_agcc = zeros(1,K);


for c=1:K
    stats=permutation_htest2_np([P2emp(:,c)',Pstatessimul(I,c)'],[ones(1,numel(P2emp(:,c))) 2*ones(1,numel(Pstatessimul(I,c)))],1000,0.05,'one_sample_ttest');
    Pval_ctrl(c)=min(stats.pvals);

    stats=permutation_htest2_np([P1emp(:,c)',Pstatessimul(I,c)'],[ones(1,numel(P1emp(:,c))) 2*ones(1,numel(Pstatessimul(I,c)))],1000,0.05,'one_sample_ttest');
    Pval_agcc(c)=min(stats.pvals);
end

%% plot PMS, simulated and empirical

figure
for c=1:K
    subplot(1,K,c)
    Group1=squeeze(Pstatessimul(I,c)); % simulated 
    Group2=squeeze(P2emp(:,c)); % empirical control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    if c==5
        title({['Simulated woCC vs. Empirical Control'] ['State #' num2str(c)]}, 'FontSize', 12);
    else
        title({['State #' num2str(c)]}, 'FontSize', 12);
    end
    
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'sim', 'emp'},'FontSize', 12)
    if Pval_ctrl(K)<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
end

figure
for c=1:K
    subplot(1,K,c)
    Group1=squeeze(Pstatessimul(I,c)); % simulated 
    Group2=squeeze(P1emp(:,c)); % empirical control
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    if c==5
        title({['Simulated woCC vs. Empirical AgCC'] ['State #' num2str(c)]}, 'FontSize', 12);
    else
        title({['State #' num2str(c)]}, 'FontSize', 12);
    end
    
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'sim', 'emp'},'FontSize', 12)
    if Pval_agcc(K)<0.05
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if c==1
        ylabel('Probability','FontSize', 12)
    end
    box off
end
