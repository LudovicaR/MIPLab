%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANALYSIS OF THE LEIDA RESULTS
%
% Written by
% Ludovica Romanin - june 2020
% 
% This script plots the results obtained after the application of LEiDA to
% the datasets. It also performs t-tests on the PMS space metrics, to
% assess group differences. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% load LEiDA results 
[file,path] = uigetfile(['../../results/leida/LEiDA_results','*.mat']);
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end
 
load(['../../results/leida/',file]);

addpath('../functions')

%% For plots, put areas removed back to zero

% load file containing the index of regions with timecourses at zero
% this is used to remove these regions from the subsequent analysis 
if size(Controls_TCS, 1) < 200
    load ../../data/TCS/areas_zero.mat;
else
    load ../../data/TCS/areas_zero_246.mat
end

V = Vemp;

% set areas that were removed back to zero for the plots
zero_vec = zeros(size(All_Subjs_TCS(1,:,:)));
for i=1:length(areas_zero)
    area = areas_zero(i);
    All_Subjs_TCS = cat(1,All_Subjs_TCS(1:area-1,:,:),zero_vec,All_Subjs_TCS(area:end,:,:));
    V = cat(2,V(:,1:area-1),zeros(1,size(V,1))',V(:,area:end));
end

N_areas = size(All_Subjs_TCS,1);

% put labels back at original number of areas
labels = 1:size(AgCC_TCS,1);

%% regions involved in the different states

% this section can only be run when using Desikan atlas 
atlas = importdata('../../data/atlas/dbs80symm_labels.txt'); 
regions_pos = {};
regions_neg = {};

for i=1:K
    regions_pos{i} = atlas(labels(find(V(i,:)>0)));
    regions_neg{i} = atlas(labels(find(V(i,:)<0)));
end

% regions at zero
regions_zero = atlas(labels(find(V(1,:)==0)));

%% test significance of complete vs. partial AgCC

% indeces of partial and complete AgCC
idx_partial = [1,2,10,11,12];
idx_complete = [3,4,5,6,7,8,9,13];

LTemp_partial=squeeze(LT1emp(idx_partial,:));
LTemp_complete=squeeze(LT1emp(idx_complete,:));

%% complete vs. partial AgCC

[sig_unc_agcc_higher_p,sig_corr_agcc_higher_p]=permutation_test([Pemp_partial',Pemp_complete']',numel(idx_partial),0.05,19,'mean','1');
[sig_unc_agcc_higher_lt,sig_corr_agcc_higher_lt]=permutation_test([LTemp_partial',LTemp_complete']',numel(idx_partial),0.05,19,'mean','1');

[sig_unc_agcc_lower_p,sig_corr_agcc_lower_p]=permutation_test([Pemp_partial',Pemp_complete']',numel(idx_partial),0.05,19,'mean','2');
[sig_unc_agcc_lower_lt,sig_corr_agcc_lower_lt]=permutation_test([LTemp_partial',LTemp_complete']',numel(idx_partial),0.05,19,'mean','2');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_partial(:,c)); % partial AgCC
    Group2=squeeze(Pemp_complete(:,c)); % complete AgCC
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'pAgCC', 'cAgCC'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
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
    set(gca,'XTickLabel',{'pAgCC', 'cAgCC'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])

clear sig_unc_agcc_higher_p sig_corr_agcc_higher_p sig_unc_agcc_higher_lt sig_corr_agcc_higher_lt;
clear sig_unc_agcc_lower_p sig_corr_agcc_lower_p sig_unc_agcc_lower_lt sig_corr_agcc_lower_lt;

%% test significance of AgCC vs. Control 

[sig_unc_agcc_higher_p,sig_corr_agcc_higher_p]=permutation_test([P1emp', P2emp']',numel(P1emp(:,1)),0.05,19,'mean','1');
[sig_unc_agcc_higher_lt,sig_corr_agcc_higher_lt]=permutation_test([LT1emp', LT2emp']',numel(LT1emp(:,1)),0.05,19,'mean','1');

[sig_unc_agcc_lower_p,sig_corr_agcc_lower_p]=permutation_test([P1emp', P2emp']',numel(P1emp(:,1)),0.05,19,'mean','2');
[sig_unc_agcc_lower_lt,sig_corr_agcc_lower_lt]=permutation_test([LT1emp', LT2emp']',numel(LT1emp(:,1)),0.05,19,'mean','2');


for c=1:K
    subplot(2,K,c)
    Group1=squeeze(P(1:NSUB_AgCC,k,ind_sort(c))); % AgCC
    Group2=squeeze(P(NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls,k,ind_sort(c))); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'AgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
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
    set(gca,'XTickLabel',{'AgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])

clear sig_unc_agcc_higher_p sig_corr_agcc_higher_p sig_unc_agcc_higher_lt sig_corr_agcc_higher_lt;
clear sig_unc_agcc_lower_p sig_corr_agcc_lower_p sig_unc_agcc_lower_lt sig_corr_agcc_lower_lt;


%% test significance of AgCC complete vs. Control

[sig_unc_agcc_higher_p,sig_corr_agcc_higher_p]=permutation_test([Pemp_complete', P2emp']',numel(Pemp_complete(:,1)),0.05,999,'mean','1');
[sig_unc_agcc_higher_lt,sig_corr_agcc_higher_lt]=permutation_test([LTemp_complete', LT2emp']',numel(LTemp_complete(:,1)),0.05,999,'mean','1');

[sig_unc_agcc_lower_p,sig_corr_agcc_lower_p]=permutation_test([Pemp_complete', P2emp']',numel(Pemp_complete(:,1)),0.05,999,'mean','2');
[sig_unc_agcc_lower_lt,sig_corr_agcc_lower_lt]=permutation_test([LTemp_complete', LT2emp']',numel(LTemp_complete(:,1)),0.05,999,'mean','2');


for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_complete(:,c)); % AgCC complete
    Group2=squeeze(P2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'cAgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
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
    set(gca,'XTickLabel',{'cAgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])

clear sig_unc_agcc_higher_p sig_corr_agcc_higher_p sig_unc_agcc_higher_lt sig_corr_agcc_higher_lt;
clear sig_unc_agcc_lower_p sig_corr_agcc_lower_p sig_unc_agcc_lower_lt sig_corr_agcc_lower_lt;

%% test significance of AgCC partial vs. Control

[sig_unc_agcc_higher_p,sig_corr_agcc_higher_p]=permutation_test([Pemp_partial', P2emp']',numel(Pemp_partial(:,1)),0.05,999,'mean','1');
[sig_unc_agcc_higher_lt,sig_corr_agcc_higher_lt]=permutation_test([LTemp_partial', LT2emp']',numel(LTemp_partial(:,1)),0.05,999,'mean','1');

[sig_unc_agcc_lower_p,sig_corr_agcc_lower_p]=permutation_test([Pemp_partial', P2emp']',numel(Pemp_partial(:,1)),0.05,999,'mean','2');
[sig_unc_agcc_lower_lt,sig_corr_agcc_lower_lt]=permutation_test([LTemp_partial', LT2emp']',numel(LTemp_partial(:,1)),0.05,999,'mean','2');

for c=1:K
    subplot(2,K,c)
    Group1=squeeze(Pemp_partial(:,c)); % AgCC complete
    Group2=squeeze(P2emp(:,c)); % Controls
    bar([mean(Group1) mean(Group2)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    title({['State #' num2str(c)]}, 'FontSize', 12)
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'pAgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_p(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_p(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
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
    set(gca,'XTickLabel',{'pAgCC', 'Ctrl'},'FontSize', 12)
    xtickangle(-60)
    if sig_unc_agcc_higher_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_higher_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if sig_unc_agcc_lower_lt(c) == 1
        plot(1.5,max([mean(Group1) mean(Group2)])+.01,'*k')
    end
    if sig_corr_agcc_lower_lt(c)==1
        plot(1.5,max([mean(Group1) mean(Group2)])+.02,'*r')
    end
    if c==1
        ylabel('Lifetime (seconds)','FontSize', 12)
    end
    box off
    
end

x0=10;
y0=10;
width=1500;
height=400;
set(gcf,'position',[x0,y0,width,height])


clear sig_unc_agcc_higher_p sig_corr_agcc_higher_p sig_unc_agcc_higher_lt sig_corr_agcc_higher_lt;
clear sig_unc_agcc_lower_p sig_corr_agcc_lower_p sig_unc_agcc_lower_lt sig_corr_agcc_lower_lt;

%% functional connectivity matrices for each states
for c=1:K
    figure()
    FC_V=V(c,1:N_areas)'*V(c,1:N_areas);
    colormap(jet)
    imagesc(FC_V)
    colorbar
    axis square
    title({['State #' num2str(c)]}, 'FontSize', 12)
    ylabel('Brain area #')
    xlabel('Brain area #')
    %saveas(gcf, ['../../figures/K18_FC_state' num2str(c) '.png'])
end

