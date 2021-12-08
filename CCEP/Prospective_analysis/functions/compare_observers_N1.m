function compare_observers_N1(myDataPath)
% function is used to compare the visual checks of the N1's of two
% observers. 
% Checked files are compared and the result it one file only containing
% N1's that were scored as N1 by BOTH observers.

files = dir(fullfile(myDataPath.CCEP_allpat));

dataBase = struct;                      

for i = 1:size(ccep_allPat.sub_labels,2)
    dataBase(i).sub_label = ccep_allPat.sub_labels{i}; 
    
     % Find rows with the sub_label of interest 
    respLoc = find(contains({files(:).name},ccep_allPat.sub_labels{i}));        
  
     % load all both the SPESclin and SPESprop of the patient
    for j=1:size(respLoc,2)                                                      
       if contains(files(respLoc(j)).name,'clin_reref_check.')               % Change to load the files of interest 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_clin = ccep_clin;
          dataBase(i).filenameClin = files(respLoc(j)).name;
            
       elseif contains(files(respLoc(j)).name,'prop_reref_check.')           % Change to load the files of interest
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep_prop = ccep_prop;   
          dataBase(i).filenameProp = files(respLoc(j)).name;
       end
    end
end



end