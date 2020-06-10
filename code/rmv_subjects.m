function rmv_subjects()
    load AgCC_TCS.mat
    
    nb_sub = [3 9 12];
    % nb_sub: list with the subjects' number to remove from the dataset
    % found from inspection after preprocessing 
    
    AgCC_TCS(:,:,nb_sub) = [];

    save AgCC_TCS_new.mat AgCC_TCS
end