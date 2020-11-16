%% Brain plots
%% Import data and specify parameters

subcortical = 1; %include subcortical areas 
saturate = 1;

load('../../results/LEiDA_results_K_16_N800.mat');

%% For plots, put areas removed back to zero
V = Vemp;

% set areas that were removed back to zero for the plots
zero_vec = zeros(size(All_Subjs_TCS(1,:,:)));
for i=1:length(areas_zero)
    area = areas_zero(i);
    All_Subjs_TCS = cat(1,All_Subjs_TCS(1:area-1,:,:),zero_vec,All_Subjs_TCS(area:end,:,:));
    V = cat(2,V(:,1:area-1),zeros(1,size(V,1))',V(:,area:end));
end

n_ROI = size(All_Subjs_TCS,1);

% put labels back at 80
labels = 1:819;

%% choose state to plot
state = 16;
CC2 = V(state,1:n_ROI);

%% Atlas + Selection

CM=zeros(n_ROI,n_ROI);

%CodeBookpath=which(strcat(parc{parcellation},'_subc_codebook.mat'));


CodeBookpath=strcat(pwd,'/Schaefer800_subc_codebook.mat');
CodeBook=load(CodeBookpath);
CodeBook=CodeBook.codeBook;

clear CodeBook2
if subcortical==0
    for i=1:size(CodeBook.id,1)-19
        CodeBook2.id(i,1)=CodeBook.id(i);
        CodeBook2.name{i,1}=CodeBook.name{i};
        CodeBook2.sname{i,1}=CodeBook.sname{i};
        CodeBook2.rname{i,1}=CodeBook.rname{i};
        CodeBook2.center{i}=CodeBook.center{i};
        CodeBook2.voxels(i)=CodeBook.voxels(i);    
    end
    CodeBook2.num=size(CodeBook.id,1)-19;
    CodeBook=CodeBook2;
end

%% adjust Cvalues for saturation (to eliminate outliers peaks)
if saturate
thr=1;
CC2new=CC2;
CC2new(find(CC2>thr))=0;
CC2new(find(CC2>thr))=max(CC2new);
CC2new(find(CC2<-thr))=0;
CC2new(find(CC2<-thr))=min(CC2new);
CC2=CC2new;
end

%% plot with normal color scheme 
CC=abs(CC2)*10;
T_conn=0.9;
Factor_SphereSize=max(CC);
Factor_Col=1;
Exp_Sphere=1;

Colormap_edges='jet';
Colormap_nodes='jet';

Gamma=0.5;
LinearWeight=0.3;
CA=[-1 1]; 


%%
View=1;
PlotBrainGraph(CM,CC,CC2,CodeBook,T_conn,Factor_SphereSize,...
    Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
    LinearWeight,CA)
saveas(gcf, strcat('../../figures/state',int2str(state),'_axial.png'));

View=2;
PlotBrainGraph(CM,CC,CC2,CodeBook,T_conn,Factor_SphereSize,...
    Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
    LinearWeight,CA)
saveas(gcf, strcat('../../figures/state',int2str(state),'_sagittal.png'));

View=3;
PlotBrainGraph(CM,CC,CC2,CodeBook,T_conn,Factor_SphereSize,...
    Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
    LinearWeight,CA)
saveas(gcf, strcat('../../figures/state',int2str(state),'_coronal.png'));

hold on
