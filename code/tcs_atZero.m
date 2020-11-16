% Script to create a .txt with the subject, area and group having a
% BOLD timecourse at zero. 
% Written on March 2020

clear all;

%load data matrices for the two groups
%load AgCC_TCS_new.mat
%load Controls_TCS.mat

% Schaefer-800
load AgCC_TCS_800.mat
load Controls_TCS_800.mat

%parameters
N_areas=size(Controls_TCS,1); % this is the number of brain regions / parcellations
TP=size(Controls_TCS,2); % volumes of the fMRIs
NSUB_Controls=size(Controls_TCS,3);
NSUB_AgCC=size(AgCC_TCS,3);
%create unique variable with 2 groups concatenated
All_Subjs_TCS=zeros(N_areas,TP,NSUB_Controls+NSUB_AgCC);
All_Subjs_TCS(:,:,1:NSUB_AgCC)=AgCC_TCS;%from 1 to 13
All_Subjs_TCS(:,:,NSUB_AgCC+1:NSUB_AgCC+NSUB_Controls)=Controls_TCS;%from 14 to 41
NSUB=size(All_Subjs_TCS,3);

% counter to have a different field for each subject-area pair
nb = 1;

% lists to display the subjects and number of areas per subjects with timecourses zero 
subjects_zero = [];
areas_zero = [];

AgCC_zero = 0;
Control_zero = 0;

for s=1:NSUB
    % subject's BOLD signals
    BOLD = All_Subjs_TCS(:,:,s);
    
    areas = []; 
    labels = [];
    
    for seed=1:N_areas
        
        labels = [labels, seed];
        
        if BOLD(seed,:) == zeros(TP)
            
            % create a structure 
            timeZero(nb).subject = s;
            
            timeZero(nb).area = seed;
            areas = [areas seed];
            
            if s <=NSUB_AgCC
                timeZero(nb).group = 'AgCC';
                AgCC_zero = AgCC_zero+1;
            else
                timeZero(nb).group = 'Control';
                Control_zero = Control_zero+1;
            end
            
            % increment the counter
            nb = nb+1;
        end
    end
    
    if length(areas) ~= 0
        subjects_zero = [subjects_zero, s];
    end
    areas_zero = [areas_zero, areas];

end

areas_zero = unique(areas_zero.','rows').';

labels(areas_zero) = [];

% save the structure in a .txt file
writetable(struct2table(timeZero), 'timecourses_zero_800.txt')
save areas_zero_800.mat areas_zero labels

