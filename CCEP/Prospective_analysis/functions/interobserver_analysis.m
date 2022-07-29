function dataBase = interobserver_analysis(myDataPath)
% Determine the Cohen's Kappa interobserver variability
% Determine this with checked files of two raters/observers
files = dir(fullfile(myDataPath.CCEP_interObVar));

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
kappa = zeros(size(uni_sublabel,2),2);  % preallocation
dataBase = struct;
for subj = 1: size(uni_sublabel,2)  
    for task= 1:size(uni_tasklabel,2)
    % kan waarschijnlijk sneller door find te gebruiken en dan binair op te
    % slaan en in te laden

    fun = @(s)~cellfun('isempty',strfind({files.name}',s));                 %#ok<STRCL1> 

    pattern_R1 = {uni_sublabel{subj}, uni_tasklabel{task}, '_SB'};
    pattern_R2 = {uni_sublabel{subj}, uni_tasklabel{task}, '_DvB'};
    
    out_R1 = cellfun(fun,pattern_R1,'UniformOutput',false);
    out_R2 = cellfun(fun,pattern_R2,'UniformOutput',false);

    idx_rater1 = all(horzcat(out_R1{:}),2);
    idx_rater2 = all(horzcat(out_R2{:}),2);

    rater1 = load(fullfile(files(idx_rater1).folder, files(idx_rater1).name));
    rater2 = load(fullfile(files(idx_rater2).folder, files(idx_rater2).name));

    disp(files(idx_rater1).name)
    disp(files(idx_rater2).name)

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
    C = confusionmat(R1_obs, R2_obs);          %disp(C);             % Convert to confusion matrix
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


    %% Determine mean N1-latency of both observers
    % A new file is created to make sure that the original observer values
    % remain.

    filefolder = fullfile(myDataPath.CCEP_interObVar);
    filename = ['sub-',uni_sublabel{subj},'_ses-1_',['task-',uni_tasklabel{task}],'_meanN1Lat.mat'];

    % When observers both approved the CCEP though selected another peak
    % more than 5 samples apart, then these responses will be plotted and
    % the correct N1-peak can be selected.
    if isequal(uni_tasklabel{task} , 'SPESclin')
        if ~exist(fullfile(filefolder,filename),'file')
            ccep_clin = select_n1_latency(rater1, rater2, myDataPath, uni_sublabel{subj},uni_tasklabel{task}, subj);

        else    % When the file _meanN1Lat.mat already exist, load that one.
            clin = load(fullfile(filefolder,filename),'dataBase');
            ccep_clin = clin.dataBase;
        end

        dataBase(subj).ccep_clin = ccep_clin.ccep_clin;
        dataBase(subj).ccep_clin.Ckappa = kappa(subj,task);
        dataBase(subj).ccep_clin.tb_channels = rater1.ccep.tb_channels; % Does not matter to take rater 1 or rater 2
        dataBase(subj).ccep_clin.n1_sample_ori = rater1.ccep.n1_peak_sample ;
        dataBase(subj).ccep_clin.n1_amp_ori = rater1.ccep.n1_peak_amplitude ;

    elseif isequal(uni_tasklabel{task} , 'SPESprop')
        
        if ~exist(fullfile(filefolder,filename),'file')
            ccep_prop = select_n1_latency(rater1, rater2, myDataPath, uni_sublabel{subj},uni_tasklabel{task}, subj);

        else
            prop = load(fullfile(filefolder,filename),'dataBase');
            ccep_prop = prop.dataBase;
        end

        dataBase(subj).ccep_prop = ccep_prop.ccep_prop;
        dataBase(subj).ccep_prop.Ckappa = kappa(subj,task);
        dataBase(subj).ccep_prop.tb_channels = rater1.ccep.tb_channels;
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

fileName = 'kappa_scores.xlsx';
writetable(T_kappa  ,[targetFolder, fileName])


%% Remove bad or silicon electrodes when they are in clin or prop
% Stimulation pairs or electrodes that met one of the exclusion criteria 
% (i.e. noisy, on top of other electrodes) were excluded in both the 
% clinical-SPES and propofol-SPES

for pat = 1:size(dataBase,2)
   if any(ismember(dataBase(pat).ccep_clin.tb_channels.status , 'bad')) || any(ismember(dataBase(pat).ccep_prop.tb_channels.status,'bad'))
        bad_clin = find(ismember(dataBase(pat).ccep_clin.tb_channels.status , 'bad'));
        bad_prop = find(ismember(dataBase(pat).ccep_prop.tb_channels.status,'bad'));
        bad_chan = unique([bad_clin;bad_prop]);

        % Remove responses of bad channels in clin and prop
        if length(dataBase(pat).ccep_clin.tb_channels.status) == length(dataBase(pat).ccep_clin.n1_peak_sample)
            dataBase(pat).ccep_clin.n1_peak_sample(bad_chan,:) = NaN;
            dataBase(pat).ccep_prop.n1_peak_sample(bad_chan,:) = NaN;
        else 
            sprintf('length of N1 matrix is not equal to the number of channels in tb_channels')
        end

        % Remove bad electrodes when they were part of stimulation pair
        % Normally bad electrodes are not stimulated, though an extra check
        % Does not matter whether you check the stimp-location of the
        % electrode in clin or in prop, the bad_chan already made sure to
        % check clin and prop.
        stimp1 = find(ismember(dataBase(pat).ccep_clin.stimsets_avg(:,1), bad_chan));
        stimp2 = find(ismember(dataBase(pat).ccep_clin.stimsets_avg(:,2), bad_chan));

        if ~isempty(stimp1) || ~isempty(stimp2) % When either one of the stimpair electrodes was a bad electrode
            bad_stimp = unique([stimp1,stimp2]);
            dataBase(pat).ccep_clin.n1_peak_sample(:,bad_stimp) = NaN;
            dataBase(pat).ccep_prop.n1_peak_sample(:,bad_stimp) = NaN;

        else

        end


   end



    %% Check for depth electrodes and remove them
    idx_depth = ismember(dataBase(pat).ccep_clin.tb_channels.group, 'depth');

    if sum(idx_depth) > 0
        dataBase(pat).ccep_clin.n1_peak_sample(idx_depth,:) = NaN;
        dataBase(pat).ccep_prop.n1_peak_sample  (idx_depth,:) = NaN;  
               
    end
end

end