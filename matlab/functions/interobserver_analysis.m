% combine the results of the two observers:
% - CCEPs are correct, if they are annotated by both observers
% - N1-latencies are corrected: when the difference between both observers
% is >5 samples: the N1-latency is annotated again. Otherwise, the average
% between both annotated N1-latencies is calculated. 

% INPUT:
% - rater1
%   struct containing the following fields:
%   - n1_peak_amplitude_check
%   matrix[channels x stim pairs] with annotated amplitudes of N1-peaks
%   - n1_peak_sample_check
%   matrix[channels x stim pairs] with annotated samples of N1-peaks
%   - n1_peak_amplitude
%   matrix[channels x stim pairs] with amplitudes of detected N1-peaks
%   - n1_peak_sample
%   matrix[channels x stim pairs] with samples of detected N1-peaks
%   - stimchans_avg
%   cell[stim pairs x 2] with names of stimulated electrodes
%   - stimpnames_avg
%   cell[1 x stim pairs] with names of stimulated electrodes, for example 'C01-C02'
%   - stimsets_avg 
%   matrix[stim pairs x 2] with numbers of stimulated electrodes
%   - dataName
%   string with filename of eeg in which SPES-events are detected and
%   annotated.

% - rater2: similar struct as rater1, but from another observer

% - myDataPath
%   struct with field
%   - CCEPpath: containing the path where to save derivatives. 

function interobserver_analysis(rater1, rater2, myDataPath)

%% Check where rater1 and rater2 have different observations

diff_obs = find((~isnan(rater1.n1_peak_sample_check) + ~isnan(rater2.n1_peak_sample_check))==1);

%% Exclude different observations from the analysis

rater1.n1_peak_sample_check(diff_obs) = NaN;  
rater1.n1_peak_amplitude_check(diff_obs) = NaN;

rater2.n1_peak_sample_check(diff_obs) = NaN;  
rater2.n1_peak_amplitude_check(diff_obs) = NaN;

% As a result, the matrices for rater1 and rater2 should be equal.
if ~isequal(~isnan(rater1.n1_peak_sample_check),~isnan(rater2.n1_peak_sample_check))
    error('Matrices of rater1 and rater2 are not equal after removal of different observations. Please check manually!')
end

%% Determine correct N1-latency when observers annotated a different N1-peak.
% When observers both scored a CCEP but the selected peaks were more than 5
% samples apart, then these responses will be plotted and the correct
% N1-peak can be selected.

filefolder = fullfile(myDataPath.CCEPpath,'checkedN1s');
[~,filename] = fileparts(rater1.dataName);
replaceRun = strfind(filename,'run');
filename = [filename(1:replaceRun-1), 'N1sChecked_comb.mat'];

% only run this if the file
% sub-PRIOSXX_ses-X_task-SPESXXXX_N1sChecked_comb.mat does not exist yet.
if ~exist(fullfile(filefolder,filename),'file')
    
    select_n1_latency(rater1, rater2, myDataPath);

end

end % end function