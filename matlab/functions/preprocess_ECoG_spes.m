% this script divides the data into epochs time-locked to the stimulus
% pulse. It takes into account: seizures, burst suppression and pulses must
% be applied with an interval of at least 3s. Each stimulus pair must be
% stimulated in anodal and cathodal direction (so C01-C02 and C02-C01),
% otherwise, it is removed from the data.

% INPUTS: 
% - cfg
%   is a struct containing the following fields:
%   - epoch_length
%   scalar defining the total length that each epoch should be
%   - epoch_prestim
%   scalar defining the length of the data pre-stimulus
%   - minstim
%   scalar defining the minimal number of stimuli that should be applied to
%   each stimulation pair.

% - dataBase
%   is a struct containing the following fields:
%   - tb_events
%   table containing information regarding events in the eeg-file (see BIDS
%   structure)
%   - ch
%   cell[channels x 1] containing the channel names 
%   - sub_label
%   string containing the name of the subject
%   - ccep_header
%   struct containing, among other things, the field: Fs (sample frequency)
%   - data
%   matrix[channels x samples] containg the un-epoched eeg-signals 

% OUTPUTS:
% - dataBase: the following fields are added in this struct
%   - cc_stimsets
%   matrix[stim pairs x 2] containing all the electrode numbers that are
%   stimulated (C1-C2 and C2-C1 are combined)
%   - cc_stimchans
%   cell[stim pairs x 2] containing all the electrode names that are
%   stimulated (C1-C2 and C2-C1 are combined)
%   - cc_epoch_sorted
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined)
%   - cc_epoch_sorted_avg
%   matrix[channels x stimulus pairs x samples] containing averaged
%   responses to stimulation (averaged from cc_epoch_sorted)
%   - tt_epoch_sorted
%   matrix[trials x stimulus pairs x samples] containing epoched time
%   points
%   - tt
%   vector[1 x samples] containing epoched time points based on the general
%   epoch length and time window pre-stimulus
%   - Burstsup
%   vector[1 x samples] containing all samples that are within periods when
%   burst suppression is observed
%   - Seizure
%   vector[1 x samples] containing all samples that are within periods when
%   a seizure is observed

function dataBase = preprocess_ECoG_spes(dataBase,cfg)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;
minstim = cfg.minstim;

for subj = 1:size(dataBase,2)
    
    %% Determine Burst suppression periods to enable removing stimuli in these periods (test in PRIOS03(?))
    BS_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.trial_type,'burst_suppression'));
    BS_stop = dataBase(subj).tb_events.sample_end(strcmp(dataBase(subj).tb_events.trial_type,'burst_suppression'));

    if iscell(BS_start)
        BS_start = str2double(BS_start);
    end
    if iscell(BS_stop)
        BS_stop = str2double(BS_stop);
    end

    BS_sample_tmp = cell(1,size(BS_start,1));          
    for iBS = 1:size(BS_start,1)
        BS_sample_tmp{1,iBS} = BS_start(iBS):BS_stop(iBS);  
    end

    BS_sample = horzcat(BS_sample_tmp{:});
    dataBase(subj).Burstsup = BS_sample;

    clear BS_sample_tmp BS_start BS_stop

    %% Determine seizure periods to enable removing stimuli in these periods (test in PRIOS01)
    % This is yet no electrode specific
    SZ_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.trial_type,'seizure'));
    SZ_stop = dataBase(subj).tb_events.sample_end(strcmp(dataBase(subj).tb_events.trial_type,'seizure'));

    if iscell(SZ_start)
        SZ_start = str2double(SZ_start);
    end
    if iscell(SZ_stop)
        SZ_stop = str2double(SZ_stop);
    end

    SZ_sample_tmp = cell(1,size(SZ_start,1));
    for iSZ = 1:size(SZ_start,1)
        SZ_sample_tmp{1,iSZ}  = SZ_start(iSZ):SZ_stop(iSZ);              
    end

    SZ_sample = horzcat(SZ_sample_tmp{:});
    dataBase(subj).Seizure = SZ_sample;

    clear SZ_sample_tmp SZ_start SZ_stop

    %% Remove stimulus pairs with less than minimal interstim time
    % Stimulations with too little interstimulus time need to be removed
    % since it cannot be guaranteed that the signal went back to baseline
    % before the next stimulus was given.

    idx_tbstim = contains(dataBase.tb_events.trial_type,'electrical_stimulation');
    tb_stim = dataBase.tb_events(idx_tbstim,:);

    % Keep stimulations that are at least 3 seconds apart (because of
    % mechanical problems with the stimulator
    idx_keep5s = [true(1); diff([tb_stim.onset])>=3];
    tb_stim = tb_stim(idx_keep5s,:);

    % keep stimulations that are not applied during a seizure
    idx_keepSz = ~ismember(tb_stim.sample_start,SZ_sample);
    tb_stim = tb_stim(idx_keepSz,:);

    % keep stimulations that are not applied during burst suppression
    idx_keepBS = ~ismember(tb_stim.sample_start,BS_sample);
    tb_stim = tb_stim(idx_keepBS,:);
    
    %% Unique stimulation pairs
    stimpair = tb_stim.electrical_stimulation_site( ...
        contains(tb_stim.sub_type,'SPES') & ...
        ~contains(tb_stim.electrical_stimulation_site,'n/a'));

    stimelek = NaN(size(stimpair,1),2);
    for stimp = 1:size(stimpair,1)
        stimchans = strsplit(stimpair{stimp},'-');
        for chan = 1:2
            stimelek(stimp,chan) = find(strcmp(stimchans{chan}, dataBase(subj).ch)==1);
        end
    end
    
    clear stimchans

    % determine stimpairs with direction (C01-C02, and C02-C01)
    [cc_stimsets_dir,~,IC_dir] = unique(stimelek,'rows');
    % determine stimpairs without direction (C01-C02 and C02-C01 = C01-C02)
    [cc_stimsets_comb,~,IC_comb] = unique(sort(stimelek,2),'rows');

    %% Remove stimulus pairs with less than minimum number of stimulations
    n = histcounts(IC_dir,'BinMethod','integers');

    if any(diff(n) ~= 0) % if any pair is stimulated a different amount

        if n < minstim
            stimremove = n<minstim;
            % remove stim pairs in both directions
            remove_elec = cc_stimsets_dir(stimremove,:);
            remove_stimp = find(cc_stimsets_dir(:,2) == remove_elec(:,2) & ...
                cc_stimsets_dir(:,1)==remove_elec(:,1) | ...
                cc_stimsets_dir(:,2)==remove_elec(:,1) &  ...
                cc_stimsets_dir(:,1)==remove_elec(:,2));

            warning('%s: stimulation pair(s) are stimulated less than all others, these are removed from further analysis\n', ...
                dataBase(subj).sub_label);

            for i = 1:length(remove_stimp)-1
                stimelek(IC_dir==remove_stimp(i),:)= NaN;
                stimelek(IC_dir==remove_stimp(i+1),:)= NaN;
            end

            remove(:,1) = isnan(stimelek(:,1));
            stimelek(remove,:) = [];

            [cc_stimsets_dir] = unique(stimelek,'rows');
        end

        warning('%s: a stimulation pair is stimulated more than others \n', ...
            dataBase(subj).sub_label)
    end

    %% Number of stimulations per stimulus pair
    % in the previous step, all stimulus pairs are found: C01-C02 and
    % C02-C01 separately. In this step, these are combined when both
    % C01-C02 and C02-C01 are stimulated. If one of both misses, then this
    % stimulus pair is removed from further analysis.

    sort_cc_stimsets_dir = sort(cc_stimsets_dir,2);
    [cc_stimsets_combdir, ~, IC_combdir] = unique(sort_cc_stimsets_dir,'rows');

    % stimpairs which are stimulated in one direction
    Ncount = histcounts(IC_combdir, size(cc_stimsets_combdir,1)); 

    if any(Ncount ~= 2) % stimpairs which are stimulated in one direction
        
        singleDir = find(Ncount ~= 2);

        remove = cell(size(singleDir,2),1);
        for n = 1:size(singleDir,2)
            remove{n} = find(IC_combdir == singleDir(n));
        end

        % Remove rows with 'single' stimpairs
        sort_cc_stimsets_dir(vertcat(remove{:}),:) = [];

        % recalculate IC_combdir
        [cc_stimsets_combdir, ~, ~] = unique(sort_cc_stimsets_dir,'rows');
    end

%% check whether any stimulus pairs are removed because it is stimulated in only one direction

if isequal(cc_stimsets_combdir,cc_stimsets_comb)

    % do nothing, because no stimulus pairs have been removed.

else

    % find which stimulus pair is removed and delete this pair
    [~, ia2] = setdiff(cc_stimsets_comb, cc_stimsets_combdir,'rows');
    
    remove = cell(size(ia2,1),1);
    for n = 1:size(ia2,1)
        remove{n} = find(IC_comb == ia2(n));
    end

    stimelek(vertcat(remove{:}),:) = [];

    % find the correct stimulus pairs again
    [cc_stimsets_comb,~,IC_comb] = unique(sort(stimelek,2),'rows');

end

    %% Determine stimulation pair names
    % pre-allocation
    cc_stimchans_comb = cell(size(cc_stimsets_comb,1),2);
    
    % determine stimulus pair names for averaged stimulus
    for stimp = 1:size(cc_stimsets_comb,1)
        for chan = 1:2
            cc_stimchans_comb{stimp,chan} = dataBase(subj).ch{cc_stimsets_comb(stimp,chan)};
        end
    end

    %% Select epochs
    
    % find maximal number of stimuli
    Ncount = histcounts(IC_comb, size(cc_stimsets_comb,1)); 
    max_stim = max(Ncount);

    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    tt = (1:epoch_length*dataBase(subj).ccep_header.Fs)/dataBase(subj).ccep_header.Fs - epoch_prestim;

    % allocation
    cc_epoch_sorted = NaN(size(dataBase(subj).data,1), size(cc_stimsets_comb,1),max_stim, t);
    tt_epoch_sorted = NaN(max_stim, size(cc_stimsets_comb,1), t);

    for elec = 1:size(dataBase(subj).data,1)                    % for all channels
        for ll = 1:size(cc_stimsets_comb,1)       % for all stimulus pairs

            eventnum = find(strcmp(tb_stim.electrical_stimulation_site, ...
                [cc_stimchans_comb{ll,1}, '-',cc_stimchans_comb{ll,2}]) | ...
                strcmp(tb_stim.electrical_stimulation_site, ...
                [cc_stimchans_comb{ll,2}, '-',cc_stimchans_comb{ll,1}]) );

            if size(eventnum,1) > max_stim
                events = max_stim;
            else
                events = size(eventnum,1);
            end

            for n = 1:events
                if tb_stim.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1< 0
                    % do nothing, because epoch starts before the start of
                    % a file (and gives an error)

                elseif ismember(tb_stim.sample_start(eventnum(n)),BS_sample)
                    % do nothing because stimulation is part of burst
                    % suppression period and therefore not thrustworthy

                elseif ismember(tb_stim.sample_start(eventnum(n)),SZ_sample)
                    % do nothing because stimulation is part of a seizure
                    % seizure periods can therefore not be scored

                else
                    cc_epoch_sorted(elec,ll,n,:) = dataBase(subj).data(elec, ...
                        tb_stim.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:tb_stim.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                    tt_epoch_sorted(n,ll,:) = tb_stim.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:tb_stim.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);

                end
            end
        end
    end

    %% Average epochs  

    cc_epoch_sorted_avg = squeeze(mean(cc_epoch_sorted,3,'omitnan'));

    %% write back to dataBase
    
    dataBase(subj).cc_stimsets = cc_stimsets_comb;
    dataBase(subj).cc_stimchans = cc_stimchans_comb;
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted;
    dataBase(subj).tt = tt;

end
end