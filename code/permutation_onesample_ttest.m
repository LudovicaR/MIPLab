%Giulia 2017
%% non parametric univariate permutation test of the difference between 1 group avaerage and 1 subject/value (non param anaolg of one sample t-test)

function [sig_meas_unc,sig_meas_corr]=permutation_onesample_ttest(data_group,data_value,alpha,nSurr,mode,c)

%%perform non parametric univariate test to assess the significance of the
%%difference between two groups (of one or more measures)

%INPUT: 
%%% - data_group: matrix subjects x measures containing the measures of the
%%% group to be tested (possibly subjects x 1, if only one measure)
%%% - data_value: matrix 1 x measures, containing the values of of
%%% reference to be compared to the group (it can be one subject values, or all zeros, if we want to test if
%%% the group mean different from 0)
%%% - alpha: significance test level
%%% - nSurr=number of surrogates, pls choose only 99,999,9999, etc.
%%% normally OK with 999
%%% - mode: type of measure to consider and test:
%%% 1)'mean',2)'central_mean':mean of the 7 central values (median+-3),3)'median'
%%% - c: contrast= '1' for one-sided group>value, '2' for one-sided group<value, '3' two-sided

%OUTPUT:
%%% - sig_meas_corr: is 1 for corrected significantly different measures, 0
%%% otherwise
%%% - (sig_meas_unc:uncorrected for multiple comp, can be added to the output if needed)

data=[data_value;data_group];
n=size(data_group,1);
nMeas=size(data_group,2);
sig_meas_unc=zeros(1,nMeas);
sig_meas_corr=zeros(1,nMeas);


for i=1:nMeas
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute real difference btw groups
    
    switch mode
    case 'mean'%1)consider mean
    diff_real(i,1)=mean(data_group(:,i))-data_value(1,i);
    case 'central_mean' %2)consider mean of central values
    m1=round(n/2);
    sdata1(:,1)=sort(data_group(:,i),'ascend');
    diff_real(i,1)=mean(sdata1(m1-3:m1+3,1))-data_value(1,i); 
    case 'median' %3)consider median
    diff_real(i,1)=median(data_group(:,i))-data_value(1,i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute surrogates by permuting subjects
    for j=1:nSurr
        r=randperm(n+1);
        data_perm(:,1)=data(r,i);
       
        switch mode
        case 'mean'%1)consider mean
        diff_surr(i,j)=mean(data_perm(1:n,1))-data_perm(n+1,1);
        case 'central_mean' %2)consider mean of central values
        sdata_perm1(:,1)=sort(data_perm(1:n,1),'ascend');
        diff_surr(i,j)=mean(sdata_perm1(m1-3:m1+3,1))-data_perm(n+1,1);
        case 'median' %3)consider median
        diff_surr(i,j)=median(data_perm(1:n,1))-data_perm(n+1,1);
        end
        
    end
diff(i,:)=[diff_surr(i,:) diff_real(i,1)];
end

if c=='2'
    diff=-diff;
end
    


%p uncorrected:
%sort differences per every measure and take alpha*N+1 max
for i=1:nMeas
    if c=='3'
        s=sort(abs(diff(i,:)),'descend');
        thr(i,1)=s((alpha/2*(nSurr+1))+1);
        if diff_real(i,1)>=thr(i,1) 
            sig_meas_unc(1,i)=1;
        elseif diff_real(i,1)<=-thr(i,1) 
            sig_meas_unc(1,i)=-1;
        end
    else
        s=sort(diff(i,:),'descend');
        thr(i,1)=s((alpha*(nSurr+1))+1);
        if diff_real(i,1)>=thr(i,1)
            sig_meas_unc(1,i)=1;
        end
    end
end


%p corrected
%find maximum per surrogate+real (within the different measures)
if c=='3'
for j=1:nSurr+1
    diff_maxima(1,j)=max(abs(diff(:,j)));
end
%sort maxima and find thr (common to every measure)
sc=sort(diff_maxima,'descend');
thrc=sc(((alpha/2)*(nSurr+1))+1);
%thrcmin=sc(size(sc,2)-((alpha/2)*(nSurr+1))+1);

for i=1:nMeas
    if diff_real(i,1)>=thrc  
        sig_meas_corr(1,i)=1;
    elseif diff_real(i,1)<=-thrc
         sig_meas_corr(1,i)=-1;
    end
end


else
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


