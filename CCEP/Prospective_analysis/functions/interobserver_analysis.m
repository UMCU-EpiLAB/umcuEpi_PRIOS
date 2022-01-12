function dataBase = interobserver_analysis(myDataPath)
% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
files = dir(fullfile(myDataPath.CCEP_interObVar));

% Extract the run_labels to match REC2Stim to PRIOS
names = {files.name};

% Find task_labels
for n = 1:size(names,2)
    files(n).subj_name = names{n}(strfind(names{n},'sub-')+4:strfind(names{n},'task-')-8);
    files(n).task_label = names{n}(strfind(names{n},'task-')+5:strfind(names{n},'task-')+12);
end

files = files(all(~cellfun(@isempty,struct2cell(files))));          % Remove all empty lines

% Find unique run_labels
uni_tasklabel = unique({files.task_label});
uni_sublabel = unique({files.subj_name});

%% Load checked ECoG files
for subj = 1: size(uni_sublabel,2)  
    for task= 1:size(uni_tasklabel,2)
    
        for i = 1:size(files)
    
            if contains(files(i).name, uni_tasklabel{task}) && contains(files(i).name, uni_sublabel{subj}) && contains(files(i).name, '_SB')  % change this to observer initials
              rater1 =  load(fullfile(files(i).folder,files(i).name)); 
        
            elseif contains(files(i).name, uni_tasklabel{task}) && contains(files(i).name, uni_sublabel{subj}) && contains(files(i).name, '_DvB')  %   % change this to observer initials
               rater2 = load(fullfile(files(i).folder,files(i).name)); 

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
               
            %% Cohens kappa unweighted is used to determine the interobserver
            % variablity
%             kappa = [];
            C = confusionmat(R1_obs, R2_obs);          disp(C);             % Convert to confusion matrix
            n = sum(C(:));                                                  % get total N
            C = C./n;                                                       % Convert confusion matrix counts to proportion of n
            r = sum(C,2);                                                   % row sum
            s = sum(C);                                                     % column sum
            expected = r*s;                                                 % expected proportion for random agree
            po = sum(diag(C));                                              % Observed proportion correct
            pe = sum(diag(expected));                                       % Proportion correct expected
            kappa(subj,task) = (po-pe)/(1-pe);                                 % Cohen's kappa
                    
            fprintf('Cohens kappa between Rater1 and Rater2 for %s during %s is k=%1.4f \n ', uni_sublabel{subj}, uni_tasklabel{task}, kappa(subj,task))
     
            %% Check where rater1 and rater 2 have different observationa
            diff_obs = find((~isnan(rater1.ccep.n1_peak_sample_check) + ~isnan(rater2.ccep.n1_peak_sample_check))==1);

            % Write to an excel to be able to further analyse when necessary
            if size(diff_obs,1) > 0
                idx_loc_diff = zeros(size(diff_obs,1),2);
                
                for i = 1:size(diff_obs,1)
                    idx_loc_diff(i,1) = ceil(diff_obs(i)/size(rater1.ccep.n1_peak_amplitude,1)); % determine stimulation pair number 
                    idx_loc_diff(i,2) = diff_obs(i) - ((idx_loc_diff(i,1)-1) * size(rater1.ccep.n1_peak_amplitude,1)); % determine electrode number
                end

                varNames = {'Stimpair','Electrode'};
                T = table(idx_loc_diff(:,1),idx_loc_diff(:,2), 'VariableNames',varNames);
                 
                targetFolder = [myDataPath.CCEP_interObVar, 'Excels_diff_obs' ,'/'];
                
                if ~exist(targetFolder,'dir')
                    mkdir(targetFolder)
                end

                fileName = ['diff_obs_',uni_sublabel{subj},'_',uni_tasklabel{task},'.xlsx'];
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
      
    
            if isequal(uni_tasklabel{task},'SPESclin')      
                dataBase(subj).ccep_clin.n1_peak_amplitude = rater1.ccep.n1_peak_amplitude_check;  % Does not matter whether rater1 or rater2 is used
                dataBase(subj).ccep_clin.n1_peak_sample = rater1.ccep.n1_peak_sample_check;        % Does not matter whether rater1 or rater2 is used
                dataBase(subj).ccep_clin.Ckappa = kappa;
                dataBase(subj).ccep_clin.sub_label = uni_sublabel{subj};
                dataBase(subj).ccep_clin.task_label= uni_tasklabel{task};
                dataBase(subj).ccep_clin.stimsets_avg = rater1.ccep.stimsets_avg;
                dataBase(subj).ccep_clin.stimpnames_avg = rater1.ccep.stimpnames_avg;
                dataBase(subj).ccep_clin.ch = rater1.ccep.ch;

            elseif isequal(uni_tasklabel{task},'SPESprop')
                dataBase(subj).ccep_prop.n1_peak_amplitude = rater1.ccep.n1_peak_amplitude_check;  % Does not matter whether rater1 or rater2 is used
                dataBase(subj).ccep_prop.n1_peak_sample = rater1.ccep.n1_peak_sample_check;        % Does not matter whether rater1 or rater2 is used
                dataBase(subj).ccep_prop.Ckappa = kappa;
                dataBase(subj).ccep_prop.sub_label = uni_sublabel{subj};
                dataBase(subj).ccep_prop.task_label= uni_tasklabel{task};
                dataBase(subj).ccep_prop.stimsets_avg = rater1.ccep.stimsets_avg;
                dataBase(subj).ccep_clin.stimpnames_avg = rater1.ccep.stimpnames_avg;
                dataBase(subj).ccep_prop.ch = rater1.ccep.ch;

            end
    
    end

end


% Save kappa in Fridge
varNames = {'SPESclin','SPESprop'};
T_kappa = table(kappa(:,1),kappa(:,2), 'VariableNames',varNames);
 
targetFolder = [myDataPath.CCEP_interObVar];

if ~exist(targetFolder,'dir')
    mkdir(targetFolder)
end

fileName = ['kappa_scores.xlsx'];
writetable(T_kappa  ,[targetFolder, fileName])

end