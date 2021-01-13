%Giulia 2013
%% non parametric univariate permutation test of the difference between 2 group avaerages (non param anaolg of t-test)

function [sig_meas_unc,sig_meas_corr]=permutation_test(data,n,alpha,nSurr,mode,c)

%%perform non parametric univariate test to assess the significance of the
%%difference between two groups (of one or more measures)

%INPUT: 
%%% - data: matrix subjects x measures containing the measures to be tested
%%% (possibly subjects x 1, if only one measure; subjects have to be ordered in group1;group2)
%%% - n: number of subjects of the first group
%%% - alpha: significance test level
%%% - nSurr=number of surrogates, es.999
%%% - mode: type of measure to consider and test:
%%% 1)'mean',2)'central_mean':mean of the 7 central values (median+-3),3)'median'
%%% - c: contrast= '1' for group1>group2, '2' for group2>group1

%OUTPUT:
%%% - sig_meas_corr: is 1 for corrected significantly different measures, 0
%%% otherwise
%%% - (sig_meas_unc:uncorrected for multiple comp, can be added to the output if needed)


nSubj=size(data,1);
nMeas=size(data,2);
sig_meas_unc=zeros(1,nMeas);
sig_meas_corr=zeros(1,nMeas);

%basing on the contrast, I invert groups, so that I will compute the
%correct difference later on
if c=='1'
    g1=1:n;
    g2=n+1:nSubj;
elseif c=='2'
    g1=n+1:nSubj;
    g2=1:n;
end


for i=1:nMeas
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute real difference btw groups
    
    switch mode
    case 'mean'%1)consider mean
    diff_real(i,1)=mean(data(g1,i))-mean(data(g2,i));
    case 'central_mean' %2)consider mean of central values
    m1=round(size(g1,2)/2);
    m2=round(size(g2,2)/2);
    sdata1(:,1)=sort(data(g1,i),'ascend');
    sdata2(:,1)=sort(data(g2,i),'ascend');
    diff_real(i,1)=mean(sdata1(m1-3:m1+3,1))-mean(sdata2(m2-3:m2+3,1)); 
    case 'median' %3)consider median
    diff_real(i,1)=median(data(g1,i))-median(data(g2,i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute surrogates by permuting subjects
    for j=1:nSurr
        r=randperm(nSubj);
        data_perm(:,1)=data(r,i);
       
        switch mode
        case 'mean'%1)consider mean
        diff_surr(i,j)=mean(data_perm(g1,1))-mean(data_perm(g2,1));
        case 'central_mean' %2)consider mean of central values
        sdata_perm1(:,1)=sort(data_perm(g1,1),'ascend');
        sdata_perm2(:,1)=sort(data_perm(g2,1),'ascend');
        diff_surr(i,j)=mean(sdata_perm1(m1-3:m1+3,1))-mean(sdata_perm2(m2-3:m2+3,1));
        case 'median' %3)consider median
        diff_surr(i,j)=median(data_perm(g1,1))-median(data_perm(g2,1));
        end
        
    end
diff(i,:)=[diff_surr(i,:) diff_real(i,1)];
end

%p uncorrected:
%sort differences per every measure and take alpha*N+1 max
for i=1:nMeas
s=sort(diff(i,:),'descend');
thr(i,1)=s((alpha*(nSurr+1))+1);
if diff_real(i,1)>thr(i,1)
    sig_meas_unc(1,i)=1;
end
end


%p corrected
%find maximum per surrogate+real (within the different measures)
for j=1:nSurr+1
    diff_maxima(1,j)=max(diff(:,j));
end
%sort maxima and find thr (common to every measure)
sc=sort(diff_maxima,'descend');
thrc=sc(((alpha)*(nSurr+1))+1);
%thrcmin=sc(size(sc,2)-((alpha/2)*(nSurr+1))+1);

for i=1:nMeas
    if diff_real(i,1)>=thrc
    sig_meas_corr(1,i)=1;
    
end
end



end


