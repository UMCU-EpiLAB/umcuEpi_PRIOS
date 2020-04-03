function dataBase = preprocess_ECoG_spes(dataBase,cfg)
epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

if exist('minstim','var') == 0
    minstim = 5;
else
    minstim = cfg.minstim;
end

for subj = 1:size(dataBase,2)   
    %% define artefact period
    
    ev_artefact_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.trial_type,'artefact'));
    ev_artefact_stop = dataBase(subj).tb_events.sample_end(strcmp(dataBase(subj).tb_events.trial_type,'artefact'));
    
    if iscell(ev_artefact_start)
        ev_artefact_start = str2double(ev_artefact_start);
    end
    if iscell(ev_artefact_stop)
        ev_artefact_stop = str2double(ev_artefact_stop);
    end
    
    ev_artefact = [];
    for i=1:size(ev_artefact_start,1)
        ev_artefact = [ev_artefact, ev_artefact_start(i):ev_artefact_stop(i)];
    end
    
    %% unique stimulation pairs
    stimpair = dataBase(subj).tb_events.electrical_stimulation_site(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a')) ;
    
    stimnum = NaN(size(stimpair,1),2);
    for stimp = 1:size(stimpair,1)
        stimchans = strsplit(stimpair{stimp},'-');
        for chan = 1:2
            stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(subj).ch)==1);
        end
    end
    
    stimcur = str2double(dataBase(subj).tb_events.electrical_stimulation_current(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a')));
    
    if strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'yes')
        stimelek = [stimnum stimcur];
    elseif strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'no')
        stimelek = stimnum;
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'yes')
        stimelek = [sort(stimnum,2) stimcur];
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'no')
        stimelek = sort(stimnum,2);
    end
    
    [cc_stimsets,~,IC] = unique(stimelek,'rows');
    
    n = histcounts(IC,'BinMethod','integers');
    
    if any(diff(n) ~= 0)
        stimremove = find(n<minstim);               % remove al stimulation pairs that are stimulated less then 5 times        
        stimelek(any(IC==stimremove,2),:) = [];
        
        [cc_stimsets,~,IC] = unique(stimelek,'rows');
        n = histcounts(IC,'BinMethod','integers');
        if any(diff(n) ~= 0)
            fprintf('ERROR: %s some stimulation pairs are stimulated less/more than all others, these are removed',dataBase(subj).sub_label)
        end
        
    end
    
    cc_stimchans = cell(size(cc_stimsets,1),2);
    
    for stimp = 1:size(cc_stimsets,1)
        for chan =1:2
            cc_stimchans{stimp,chan} = dataBase(subj).ch{cc_stimsets(stimp,chan)};
        end
        
    end
    
    max_stim = median(n);
    
    dataBase(subj).cc_stimsets = cc_stimsets;
    dataBase(subj).cc_stimchans = cc_stimchans;
    dataBase(subj).max_stim = max_stim;
    
    stimdif = find(n ~= max_stim);
    for stimp =1:size(stimdif,2)
        [cc_stimchans{stimdif(stimp),1} '-' cc_stimchans{stimdif(stimp),2}]
    end
    
    %% select epochs
    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    
    % allocation
    cc_epoch_sorted = NaN(size(dataBase(subj).data,1),dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets,1),t);
    tt_epoch_sorted = NaN(dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets,1),t); % samplenumbers for each epoch
    
    for elec = 1:size(dataBase(subj).data,1) % for all channels
        for ll = 1:size(dataBase(subj).cc_stimsets,1) % for all epochs with >4 stimuli
            if strcmp(cfg.dir,'no')
                % Find the stimulationnumbers on which a stimulation pair is stimulated. 
                eventnum1 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}])); % Positive stimulation
                eventnum2 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,2}, '-',dataBase(subj).cc_stimchans{ll,1}])); % Negative stimulation
                eventnum = [eventnum1;eventnum2];
            elseif strcmp(cfg.dir,'yes')
                eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
            end
            
            if size(eventnum,1) >dataBase(subj).max_stim
                events = dataBase(subj).max_stim;
            else
                events = size(eventnum,1);
            end
            
            for n=1:events
                
                if dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1< 0
                    % do nothing, (samplestartnumber - Fs)+1 <1 means that????
                elseif ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),ev_artefact)
                    % do nothing, because part of artefact
                else
                    cc_epoch_sorted(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                    tt_epoch_sorted(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);
                end
            end
        end
    end
    
    cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
    
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    
    fprintf('...%s has been epoched and averaged... \n',dataBase(subj).sub_label)
    
end