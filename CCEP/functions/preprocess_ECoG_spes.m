
function dataBase = preprocess_ECoG_spes(dataBase,cfg,avg_stim)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

if ~any(contains(fieldnames(cfg),'minstim')) 
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
        ev_artefact = [ev_artefact, ev_artefact_start(i):ev_artefact_stop(i)]; %#ok<AGROW>
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
%         sort_stimelek = sort(stimelek,2);
%         [cc_stimsets,~,IC] = unique(sort_stimelek,'rows');  
%         
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'yes')
        stimelek = [sort(stimnum,2) stimcur];
        
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'no')
        stimelek = sort(stimnum,2);

    end
    
%     eerst sorteren!!!
%     sort_cc_stimsets = sort(cc_stimsets,2);
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
    
    % pre-allocation
    cc_stimchans = cell(size(cc_stimsets,1),2);
    stimpnames = cell(1,size(cc_stimsets,1));
    for stimp = 1:size(cc_stimsets,1)
        for chan =1:2
            cc_stimchans{stimp,chan} = dataBase(subj).ch{cc_stimsets(stimp,chan)};
        end
        stimpnames{stimp} = [cc_stimchans{stimp,1} '-' cc_stimchans{stimp,2}]; 
        
    end
    
    max_stim = median(n);
    
    dataBase(subj).cc_stimsets = cc_stimsets;
    dataBase(subj).cc_stimchans = cc_stimchans;
    dataBase(subj).stimpnames = stimpnames;
    dataBase(subj).max_stim = max_stim;
    
    stimdif = find(n ~= max_stim);
    for stimp =1:size(stimdif,2)
        [cc_stimchans{stimdif(stimp),1} '-' cc_stimchans{stimdif(stimp),2}] %#ok<NOPRT>
    end
    
    %% select epochs
    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    tt = (1:epoch_length*dataBase(subj).ccep_header.Fs)/dataBase(subj).ccep_header.Fs - epoch_prestim;

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
            

            for n = 1:events  
                if dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1< 0
                    % do nothing, because epoch starts before the start of
                    % a file (and gives an error)

                elseif ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),ev_artefact)
                    % do nothing, because part of artefact
                else
                    cc_epoch_sorted(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                    tt_epoch_sorted(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);
                end
            end
        end
    end   
    
    %% average epochs
    
    if strcmp(cfg.dir ,'yes') 
        if any(contains(fieldnames(cfg),'dir_avg'))
            if strcmp(cfg.dir_avg,'yes')
                sort_cc_stimsets = sort(cc_stimsets,2);
                [cc_stimsets_dir,~,IC] = unique(sort_cc_stimsets,'rows');
                
                % preallocation
                cc_epoch_sorted_avg = NaN(size(cc_epoch_sorted,1),size(cc_stimsets_dir,1),size(cc_epoch_sorted,4)); % [ channels x stim pairs_both x samples]
                
                for ll = 1:max(IC)
                
                    if sum(IC == ll) == 2 % both directions are stimulated
                        
                        if ~exist('avg_stim','var') % average all events in one stimulation pair
                            cc_epoch_sorted_avg(:,ll,:) = squeeze(nanmean(squeeze(nanmean(cc_epoch_sorted(:,:,IC==ll,:),2)),2));
                        elseif isempty(avg_stim)
                            cc_epoch_sorted_avg(:,ll,:) = squeeze(nanmean(squeeze(nanmean(cc_epoch_sorted(:,:,IC==ll,:),2)),2));
                        else % average only 1 or a few events in one stimulation pair
                            cc_epoch_sorted_avg(:,ll,:) =  squeeze(nanmean(squeeze(nanmean(cc_epoch_sorted(:,avg_stim,IC==ll,:),2)),2));
                        end
                    end
                end
            elseif strcmp(cfg.dir_avg,'no')
                %[cc_stimsets_dir,~,IC] = unique(cc_stimsets,'rows'); % Hier ben ik niet zeker van (10-7-2020)
                % preallocation
                cc_epoch_sorted_avg = NaN(size(cc_epoch_sorted,1),size(cc_stimsets,1),size(cc_epoch_sorted,4)); % [ channels x stim pairs_dir x samples]
                
                if ~exist('avg_stim','var') % average all events in one stimulation pair
                    cc_epoch_sorted_avg(:,ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
                elseif isempty(avg_stim)
                    cc_epoch_sorted_avg(:,ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
                else % average only 1 or a few events in one stimulation pair
                    cc_epoch_sorted_avg(:,ll,:) = squeeze(nanmean(cc_epoch_sorted(:,avg_stim,:,:),2));
                end
            end

        else % if cfg.dir_avg does not exist
            % preallocation
            cc_epoch_sorted_avg = NaN(size(cc_epoch_sorted,1),size(cc_stimsets,1),size(cc_epoch_sorted,4)); % [ channels x stim pairs_dir x samples]
            
            if ~exist('avg_stim','var') % average all events in one stimulation pair
                cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
            elseif isempty(avg_stim)
                cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
            else % average only 1 or a few events in one stimulation pair
                cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted(:,avg_stim,:,:),2));
            end
        end
            
    elseif strcmp(cfg.dir,'no')
        if ~exist('avg_stim','var') % average all events in one stimulation pair
            cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
        elseif isempty(avg_stim)
            cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted,2));
        else % average only 1 or a few events in one stimulation pair
            cc_epoch_sorted_avg(:,1:ll,:) = squeeze(nanmean(cc_epoch_sorted(:,avg_stim,:,:),2));
        end
    end
        
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted;
    dataBase(subj).tt = tt;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    dataBase(subj).stimpnames = stimpnames;
    dataBase(subj).cc_stimsets_avg = cc_stimsets_dir;
           
    fprintf('...%s has been epoched and averaged... \n',dataBase(subj).sub_label)
    

end