function interobserverKappa(myDataPath)
clc;
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
for run = 1: size(uni_runlabel,2)
    for i = 1:size(files)
   
        % PRIOS study is rater2, REC2Stim is rater1.
        if contains(files(i).name, uni_runlabel{run}) && contains(files(i).name, 'PRIOS') 
           rater2(run,:) = load(fullfile(files(i).folder,files(i).name)) ;
        elseif contains(files(i).name, uni_runlabel{run}) && contains(files(i).name, 'REC2Stim')
           rater1(run,:) =load(fullfile(files(i).folder,files(i).name)) ;
 
        end
    end
    
    N1_peak_R1 = rater1(run).ccep.checked  ;                        % Rater 1 already transformed to 1/0
    N1_peak_R2 = rater2(run).ccep_clin.n1_peak_sample_check;
    
    % Determine whether the same stimulation pairs are used. 
    % Stimulation pairs that are 'extra' are not removed from the original
    % detection matrix in the CCEP struct. 
    if size(N1_peak_R1,2) > size(N1_peak_R2,2)
       extra_rater1 = find(~ismember(rater1(run).ccep.cc_stimsets(:,1), rater2(run).ccep_clin.stimsets_avg(:,1))) ;
       N1_peak_R1(:,extra_rater1) = [];  
       rater1(run).ccep.cc_stimsets(extra_rater1,:) = [];
    elseif size(N1_peak_R2,2) > size(N1_peak_R1,2)
        extra_rater2 =  find(~ismember(rater2(run).ccep_clin.stimsets_avg(:,1), rater1(run).ccep.cc_stimsets(:,1))) ;
        N1_peak_R2(:,extra_rater2) = []; 
        rater2(run).ccep_clin.stimsets_avg(extra_rater2,:) = [];
    end
    
    % Check whether the stimpairs are equal in both ratings. 
    if ~isequal(rater1(run).ccep.cc_stimsets(:,:) , rater2(run).ccep_clin.stimsets_avg(:,:))
        diff_stimp = find(rater1(run).ccep.cc_stimsets ~= rater2(run).ccep_clin.stimsets_avg);
        rater1(run).ccep.cc_stimsets(diff_stimp(1,:),:) = [];
        N1_peak_R1(:,diff_stimp(1,:)) = [];  

        rater2(run).ccep_clin.stimsets_avg(diff_stimp(1,:),:) = [];
        N1_peak_R2(:,diff_stimp(1,:)) = []; 

        
        % Check again
        if ~isequal(rater1(run).ccep.cc_stimsets(:,:) , rater2(run).ccep_clin.stimsets_avg(:,:))
            warning('Check the stimulation pairs since they are still not equal')
        end
  end

    % Convert matrix with sample numbers to binary matrix
    N1_peak_R2(~isnan(N1_peak_R2)) = 1;     % all non-NaNs are amplitudes, so N1s --> 1
    N1_peak_R2(isnan(N1_peak_R2)) = 0;      % all NaNs are no N1s --> 0
    
    % Convert matrix to array
    N1_peak_R1 = N1_peak_R1(:);
    N1_peak_R2 = N1_peak_R2(:);
    
    
    % Cohens kappa unweighted is used to determine the interobserver
    % variablity
    C = confusionmat(N1_peak_R1, N1_peak_R2)                        % Convert to confusion matrix
    n = sum(C(:));                                                  % get total N
    C = C./n;                                                       % Convert confusion matrix counts to proportion of n
    r = sum(C,2);                                                   % row sum
    s = sum(C);                                                     % column sum
    expected = r*s;                                                 % expected proportion for random agree
    po = sum(diag(C));                                              % Observed proportion correct
    pe = sum(diag(expected));                                       % Proportion correct expected
    kappa = (po-pe)/(1-pe);                                  % Cohen's kappa
  
    REC2Stim_label = rater1(run).ccep.dataName(strfind(rater1(run).ccep.dataName,'sub-')+4 : strfind(rater1(run).ccep.dataName,'ses-')-2);
    PRIOS_label = rater2(run).ccep_clin.dataName(strfind(rater2(run).ccep_clin.dataName,'sub-') +4 :strfind(rater2(run).ccep_clin.dataName,'ses-')-2);
    
    fprintf('Cohens kappa between Rater1 and Rater2 for %s/%s is k=%1.4f \n ', REC2Stim_label, PRIOS_label, kappa)
    
    % Find where dorien had 0 and i have 1, determine column nummer en row
    % nummer. Vind zo de elektrode en stimpair en die dan plotten in
    % pipeline01. 
    diff_rating = num2cell(find(N1_peak_R1 ~= N1_peak_R2)) ;             % find where R1 and R2 are diffferent
    num_elec = size(rater1(run).ccep.ch  ,1);                       % New stimpair after the number of electrodes
    
    stimp = 1;                                          % Stimpair_label
    for i = 1:size(diff_rating,1)
       if  (diff_rating{i}/num_elec)>1                            % For stimpairs >1
           stimp = ceil(diff_rating{i}/num_elec);
          
           diff_rating(i,2) = rater1(run).ccep.ch(diff_rating{i}-((stimp-1)*num_elec));
           diff_rating(i,3:4) = num2cell(rater2(run).ccep_clin.stimsets_avg(stimp,1:2));

       else % for stimpair 1
           diff_rating(i,2) = rater1(run).ccep.ch(diff_rating{i});
           diff_rating(i,3:4) = num2cell(rater2(run).ccep_clin.stimsets_avg(stimp,1:2));
       end
    end
    
    varNames = {'nummer','Channel','StimElec1','StimpElec2'};
    T = table(diff_rating(:,1),diff_rating(:,2),diff_rating(:,3),diff_rating(:,4), 'VariableNames',varNames);
     
    targetFolder = [myDataPath.CCEP_interObVar];
    fileName = ['Different_ratings_',PRIOS_label,'_',REC2Stim_label,'.xlsx'];
    writetable(T  ,[targetFolder, fileName])

end
end