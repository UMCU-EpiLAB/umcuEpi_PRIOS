function dataBase = interobserverKappa(dataBase, myDataPath)
% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
files = dir(fullfile(myDataPath.CCEP_interObVar));

% Extract the run_labels to match REC2Stim to PRIOS
names = {files.name};

% Find run_labels
for n = 1:size(names,2)
    files(n).run_label = names{n}(strfind(names{n},'run-'):strfind(names{n},'run-')+9);
end

files = files(all(~cellfun(@isempty,struct2cell(files))));          % Remove all empty lines

% Find unique run_labels
uni_runlabel = unique({files.run_label});

%% Load checked ECoG files
  
for runs= 1:size(uni_runlabel,2)

    for i = 1:size(files)

        if contains(files(i).name, uni_runlabel{runs}) && contains(files(i).name, '_01')  % change this to observer initials
          rater1 =  load(fullfile(files(i).folder,files(i).name)); 
    
        elseif contains(files(i).name, uni_runlabel{runs}) && contains(files(i).name, '_02')   % change this to observer initials
           rater2 = load(fullfile(files(i).folder,files(i).name)) ; 
        end
    end
        
    
    %% Convert sample and amplitude value to binary output.
        rater1.ccep.check_binary = ~isnan(rater1.ccep.n1_peak_sample_check);    % cells with a value are converted to a 1, NaNs are converted to 0      
        rater2.ccep.check_binary = ~isnan(rater2.ccep.n1_peak_sample_check);
        
        % Only the automatically detected ERs are checked
        % Find the automatically detected ERs
        if isequal(size(rater1.ccep.n1_peak_sample_check), size(rater2.ccep.n1_peak_sample_check))  % extra check whether the same number of responses is compared
            auto_det_ERs = find( ~isnan(rater1.ccep.n1_peak_sample) == 1);      % use the non-checked file! determine the difference between the non-checked file and the checked file
        else
            error('Matrix with automatically detected ERs is not equal')
        end
        
        % Determine where the detector detected an N1 since these are visually
        % checked. All signals without automatically detected N1's are not
        % visually checked and are therefore not part of the inter observer
        % agreement.
        R1_obs = rater1.ccep.check_binary(auto_det_ERs);
        R2_obs = rater2.ccep.check_binary(auto_det_ERs);
           
        % Cohens kappa unweighted is used to determine the interobserver
        % variablity
        C = confusionmat(R1_obs, R2_obs);    %disp(C);            % Convert to confusion matrix
        n = sum(C(:));                                                  % get total N
        C = C./n;                                                       % Convert confusion matrix counts to proportion of n
        r = sum(C,2);                                                   % row sum
        s = sum(C);                                                     % column sum
        expected = r*s;                                                 % expected proportion for random agree
        po = sum(diag(C));                                              % Observed proportion correct
        pe = sum(diag(expected));                                       % Proportion correct expected
        kappa = (po-pe)/(1-pe);                                         % Cohen's kappa
        
        PRIOS_label = rater2.ccep.dataName(strfind(rater2.ccep.dataName,'sub-') :strfind(rater2.ccep.dataName,'ses-')-2);
        SPES_label = rater2.ccep.dataName(strfind(rater2.ccep.dataName,'task-') +5 :strfind(rater2.ccep.dataName,'run-')-2);
        
        fprintf('Cohens kappa between Rater1 and Rater2 for %s during %s is k=%1.4f \n ', PRIOS_label, SPES_label, kappa)
 
        %% Check where rater1 and rater 2 have different observationa
        diff_obs = find((~isnan(rater1.ccep.n1_peak_sample_check) + ~isnan(rater2.ccep.n1_peak_sample_check))==1);
        
        % Write to an excel to be able to further analyse when necessary
        for i = 1: size(diff_obs,1)
            stimp_diff = ceil(diff_obs(i)/size(rater1.ccep.ch  ,1));
            diff_rat{i,1} = rater1.ccep.stimpnames_avg{stimp_diff};
            diff_rat{i,2} = rater1.ccep.ch{diff_obs(i)-((stimp_diff-1)*size(rater1.ccep.ch  ,1))};
        end
             
        if exist('diff_rat','var')              % only save when table exist
            varNames = {'Stimpair','Electrode'};
            T = table(diff_rat(:,1),diff_rat(:,2), 'VariableNames',varNames);
             
            targetFolder = [myDataPath.CCEP_interObVar];
            fileName = ['Different_ratings_clinical_',PRIOS_label,'_',SPES_label,'.xlsx'];
            writetable(T  ,[targetFolder, fileName])
        end
        
        %% Exclude different observations form the analysis
        % Clinical SPES
        rater1.ccep.n1_peak_sample_check(diff_obs) = NaN;  % Replace the values with NaN when the observer did not have the same rating
        rater2.ccep.n1_peak_sample_check(diff_obs) = NaN;

        rater1.ccep.n1_peak_amplitude_check(diff_obs) = NaN;  % Replace the values with NaN when the observer did not have the same rating
        rater2.ccep.n1_peak_amplitude_check(diff_obs) = NaN;

        % As a result, the table for Rater 1 and rater 2 should be equal. And
        % analysis can proceed with either Rater 1 or rater2 
        if ~isequal(isnan(rater1.ccep.n1_peak_sample_check), isnan(rater2.ccep.n1_peak_sample_check))
            error('Binary matrices of rater1 and rater 2 are not equal after removal of different observations.')
        end

        
        %% Save ccep's to continue with the analysis
        % I want to save the new N1-peak_check files in the dataBase file
        % to save and use for further analysis. 

%         PAT_label = rater1.ccep.dataName(strfind(rater1.ccep.dataName,'/ieeg/sub-') +6:strfind(rater1.ccep.dataName,'_ses-1')-1);% Determine patient label of current files
        loc_pat_in_dataBase = find(ismember({dataBase.sub_label}.' , PRIOS_label));       % Find where in dataBase the same pat_label is used


        if isequal(SPES_label,'SPESclin')%  SPES_label is based on the run_label
            dataBase(loc_pat_in_dataBase).ccep_clin.n1_peak_amplitude_check = rater1.ccep.n1_peak_amplitude_check;  % Does not matter whether rater1 or rater2 is used
            dataBase(loc_pat_in_dataBase).ccep_clin.n1_peak_sample_check = rater1.ccep.n1_peak_sample_check;        % Does not matter whether rater1 or rater2 is used

        elseif isequal(SPES_label,' SPESprop')
            dataBase(loc_pat_in_dataBase).ccep_prop.n1_peak_amplitude_check = rater1.ccep.n1_peak_amplitude_check;  % Does not matter whether rater1 or rater2 is used
            dataBase(loc_pat_in_dataBase).ccep_prop.n1_peak_sample_check = rater1.ccep.n1_peak_sample_check;        % Does not matter whether rater1 or rater2 is used

        end


end
end