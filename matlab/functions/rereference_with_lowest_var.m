% The trials after each stimulus are re-referenced. First, the variance is
% calculated for all signals except the stimulus pair and bad channels in
% time windows pre (-0.5:-0.1) and post stimulation (0.01:0.01). Then, the
% 10% channels with the lowest variance pre- and post stimulation are
% calculated.

% Lowest variance is used to avoid introducing CCEP's in signals that
% originally do not have a CCEP.

% INPUT:
% - dataBase:
% is a struct containing the following fields:
%   - tt
%   matrix[trials x stimulus pairs x samples] containing epoched time
%   points
%   - cc_epoch_sorted
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined)
%   - tb_channels
%   table containing information regarding the recording channels (see BIDS
%   structure)

% OUTPUTS:
% - dataBase:
% the following fields are added to this struct:
% - cc_epoch_sorted_reref
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined) --> re-referenced
% - cc_epoch_sorted_avg_reref
%   matrix[channels x stimulus pairs x samples] containing averaged
%   responses to stimulation (averaged from cc_epoch_sorted) -->
%   re-referenced
% - ref
%   matrix[stimulus pairs x trials x samples] containing the reference
%   signal that is used to re-reference the original data in
%   cc_epoch_sorted

function dataBase = rereference_with_lowest_var(dataBase)

period_preStim = find(dataBase.tt >=-0.5 & dataBase.tt <=-0.1);  % [-0.5: -0.1] pre-stim period
period_postStim = dataBase.tt >0.01 & dataBase.tt <=0.1;   % [0.01: 0.1] post-stim period
data_all = dataBase.cc_epoch_sorted;

% pre-allocation
ref_all = NaN(size(dataBase.cc_epoch_sorted,2),size(dataBase.cc_epoch_sorted,3),size(dataBase.cc_epoch_sorted,4));
cc_epoch_sorted_reref = NaN(size(dataBase.cc_epoch_sorted));

for stimp = 1:size(dataBase.cc_epoch_sorted,2)                  % For each stimulation pair.

    % Exclude bad channels and stimulated channels by making it NaN's
    bad_and_stimp = [find(strcmp(dataBase.tb_channels.status_description, ...
        'noisy (visual assessment)')); dataBase.cc_stimsets(stimp,:)'];
    data_all(bad_and_stimp,stimp,:,:) = NaN;
    
    for numstim = 1:size(data_all,3)  % for each of the trials of one stimulus pair

        % Determine median of all signals except BAD and
        % stimulated channels (common average reference: CAR)
        raw_data_stim = squeeze(data_all(:,stimp,numstim,:));
        CA_raw_data = median(raw_data_stim,'omitnan');

        % During propofol SPES there are often varying number of
        % stims per stimulation pair, and when there are only
        % NaNs in raw_data_stim, CAR_raw_data also contains
        % only NaNs. This has to be ignored in the analysis.
        if ~isnan(sum(CA_raw_data))

            % Determine variance of CAR signal PRE-stim
            var_preStimCAR = var(CA_raw_data(:, period_preStim),'omitnan');

            % Determine for each separate signal if variance
            % POST-stim is lower than variance of CAR signal
            % var(x,0 or 1) for normalization (Normalization is a
            % good technique to use when you do not know the distribution of your data)
            % var(x,1,1) = var for columns, var(x,1,2) = var for rows
            var_postStim_signal = var(squeeze(data_all(:,stimp,numstim, period_postStim)),1,2,'omitnan');

            % Determine for each separate signal if variance
            % pre-stim is lower than variance of CAR signal
            var_preStim_signal = var(squeeze(data_all(:,stimp,numstim, period_preStim )),1,2,'omitnan');

            % Keep the signals that have a variance post- AND
            % pre-stim that is lower than the variance of the
            % CAR signal Determine the median --> this is your
            % reference signal
            post_lower = find(var_postStim_signal < var_preStimCAR);
            pre_lower = find(var_preStim_signal < var_preStimCAR);
            ref_keep = post_lower(ismember(post_lower,pre_lower));

            % If there are less than 10% of the signals present
            % in the reference, Then add signals with lowest
            % variance POST-stim.
            if size(ref_keep,1) < ceil(0.1*size(data_all,1))

                [~,idx_var_post] = sort(var_postStim_signal,'ascend');                             % NAN data receives the highest scoring, and are therefore at the bottom of the ranking
                nmb_extra_signals = round(0.1*size(raw_data_stim,1))- size(ref_keep,1);            % calculate the number of extra signals to obtain 5% of the signals in the reference

                % If the lowest 10% variance is already in the ref_keep, then take the number of extra required
                % signals + the number signals that were already present in the ref_keep
                if any(ismember(ref_keep, idx_var_post(1:(size(ref_keep,1)+nmb_extra_signals))))
                    nmb_extra_signals = nmb_extra_signals + sum(ismember(ref_keep, idx_var_post(1:(size(ref_keep,1)+nmb_extra_signals))));
                end

                ref_keep2 = unique([ref_keep; idx_var_post(1:nmb_extra_signals)]);
                ref_keep = ref_keep2;
            end

            % calculate reference signal
            ref = median(squeeze(data_all(ref_keep, stimp,numstim, :)));
            ref_all(stimp,numstim,:) = ref;

            % Re-reference the individual trials. Make sure to take the original signal to make
            % sure to also rereference the bad channels and stimulated channels
            cc_epoch_sorted_reref(:,stimp,numstim,:) = squeeze(dataBase.cc_epoch_sorted(:,stimp,numstim,:)) - ref;

        else
            % do nothing because this stimulation pair is not stimulated
            % for this numstim (contains only NaN)
            cc_epoch_sorted_reref(:,stimp,numstim,:) = squeeze(dataBase.cc_epoch_sorted(:,stimp,numstim,:));
        end

    end % for-loop for each trial
end % for-loop for each stimulus pair

% Take mean of all numstims of the re-referenced signals per stimulation pair
cc_epoch_sorted_reref_avg = squeeze(mean(cc_epoch_sorted_reref,3,'omitnan'));

% write back to dataBase
dataBase.ref = ref_all;
dataBase.cc_epoch_sorted_reref = cc_epoch_sorted_reref;
dataBase.cc_epoch_sorted_reref_avg = cc_epoch_sorted_reref_avg;

fprintf('%s %s has been re-referenced. \n',dataBase.sub_label, dataBase.task_label)

end