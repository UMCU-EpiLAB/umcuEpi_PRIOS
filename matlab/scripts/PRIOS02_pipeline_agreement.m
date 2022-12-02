% This script compares the visual annotations between the two observers,
% calculates Cohen's kappa, and removes all responses that were not
% annotated by both observers. The final responses are saved and used in
% further analyses. 

clear; 
clc;

%% set paths

% add current path from folder which contains this script
rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = PRIOS_setLocalDataPath(1);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% select all subjects available (in whom 2 observers visually annotated N1-peaks)

files = dir(fullfile(myDataPath.CCEPpath,'checkedN1s'));
files = files([files(:).bytes] > 0);

if size(files,1) == 0
    error('You should first run PRIOS01_pipeline_preprocess.m by two observers before you can run this script.')
end

names = {files.name};

% Find unique sub_labels and task_labels for which you can run the consecutive codes
for n = 1:size(names,2)
    files(n).subj_name = names{n}(strfind(names{n},'sub-')+4:strfind(names{n},'task-')-8);
    files(n).task_label = names{n}(strfind(names{n},'task-')+5:strfind(names{n},'task-')+12);
end

% Find unique run_labels
all_tasklabel = unique({files.task_label});
all_sublabel = unique({files.subj_name});

% housekeeping
clear n

%% calculate Cohen's kappa

% pre-allocation
all_kappa = NaN(size(all_sublabel,2),size(all_tasklabel,2));

for subj = 1:size(all_sublabel,2)  
    for task = 1:size(all_tasklabel,2)

        %% load ECoG files

        idx_rater1 = contains(names,all_sublabel{subj}) & contains(names,all_tasklabel{task}) & contains(names,'_SB');
        idx_rater2 = contains(names,all_sublabel{subj}) & contains(names,all_tasklabel{task}) & contains(names,'_DvB');

        load(fullfile(files(idx_rater1).folder, files(idx_rater1).name));
        rater1 = ccep; % FIXTHIS: als ik de struct anders opsla, dan kan deze stap eruit en kan het direct in rater1 ingeladen worden.
        load(fullfile(files(idx_rater2).folder, files(idx_rater2).name));
        rater2 = ccep; % FIXTHIS: als ik de struct anders opsla, dan kan deze stap eruit en kan het direct in rater1 ingeladen worden.

        disp(files(idx_rater1).name)
        disp(files(idx_rater2).name)

        %% Convert sample and amplitude value to binary output.
        rater1.check_binary = ~isnan(rater1.n1_peak_sample_check);    % cells with a value are converted to a 1, NaNs are converted to 0
        rater2.check_binary = ~isnan(rater2.n1_peak_sample_check);

        % Only the automatically detected CCEPs are checked visually, so
        % find the automatically detected CCEPs
        if isequal(size(rater1.n1_peak_sample), size(rater2.n1_peak_sample))  % extra check whether the same number of responses is compared
            
            if isequal(~isnan(rater1.n1_peak_sample),~isnan(rater2.n1_peak_sample))
                auto_det_CCEPs = find(~isnan(rater1.n1_peak_sample) == 1);      % use the file with automatically detected CCEPs
            else
                error('Both observers seem to have used a different detection algorithm. Comparing both visual ratings is not possible, since only detected CCEPs are visually checked.')
            end

        else
            error('The sizes of matrices with automatically detected CCEPs differ between the two raters')
        end

        % Determine when the detector detected an N1 since only these are
        % visually checked. All signals without automatically detected N1s
        % are not visually checked and are therefore not part of the
        % inter-observer agreement.
        R1_obs = rater1.check_binary(auto_det_CCEPs);
        R2_obs = rater2.check_binary(auto_det_CCEPs);

        %% Cohens kappa unweighted 
        % is used to determine the inter-observer variablity

        kappa = cohensKappa(R1_obs, R2_obs);
        all_kappa(subj,task) = kappa;

        fprintf('Cohens unweighted kappa between Rater1 and Rater2 for %s during %s is k = %1.4f \n ', ...
            all_sublabel{subj}, all_tasklabel{task}, all_kappa(subj,task))

        %% Combine the visually annotated CCEPs of the two observers
        % Exclude responses that were scored differently between observers
        % and correct the N1-peak if observers annotated a peak with >5
        % samples difference.

        interobserver_analysis(rater1, rater2, myDataPath)

    end % for-loop tasks
end % for-loop subjects




%% moet weg
% FIXTHIS: dit moet nog weg
% %% Plot electrodes position on brain with N1-latency information
% % Determine N1-latency per brain part
% close all
% 
% plot_electrodes_on_MRI(myDataPath, dataBase)
% 
% %% Determine the distance between electrodes to determine correlation between N1-latency and electrode distance
% close all
% distance_elec_stimp(dataBase, myDataPath, av_lat_elec)
