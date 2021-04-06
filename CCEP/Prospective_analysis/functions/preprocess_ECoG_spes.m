function dataBase = preprocess_ECoG_spes(dataBase,cfg)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

if ~any(contains(fieldnames(cfg),'minstim'))
    minstim = 5;
else
    minstim = cfg.minstim;
end

% if cfg.dir_avg is not defined in config_CCEP.m, default is 'no'
if sum(contains(fieldnames(cfg),'dir_avg'))<1
    cfg.dir_avg = 'no';
end

for subj = 1:size(dataBase,2)
    
% Define artefact period - stimuli during an artefact become NaNs
% This must be defined for all electrodes or for specified electrodes.

    elec_art = find(strcmp(dataBase(subj).tb_events.trial_type,'artefact'));                % Find rows with artefact in them
    elec_art_all = find(strcmp(dataBase(subj).tb_events.electrodes_involved_onset,'all'));  % Find the artefacts which concern all electrodes
    elec_art_in_all = find(ismember(elec_art, elec_art_all));                               % done because notes like STIM_on also concern 'all'.

    % Determine the start- and stop-moment of the artefact concerning all
    % electrodes
    ev_artefact_start = dataBase(subj).tb_events.sample_start(elec_art_in_all,:);
    ev_artefact_stop = dataBase(subj).tb_events.sample_end(elec_art_in_all,:);

    if iscell(ev_artefact_start)
        ev_artefact_start = str2double(ev_artefact_start);
    end

    if iscell(ev_artefact_stop)
        ev_artefact_stop = str2double(ev_artefact_stop);
    end

    ev_artefact_all = zeros(1,100);                                                         % Preallocation with a guess of the size, structs with empty cells are removed.
    % When there are multiple artefacts, then the times samples are concatenaded
    % in the same array. So ev_artefact is not overwritten
    for i=1:size(ev_artefact_start,1)
        ev_artefact_all = [ev_artefact_all, ev_artefact_start(i):ev_artefact_stop(i)];    %#ok<AGROW>
    end

    dataBase(subj).ev_artefact_all = ev_artefact_all;
    


%% When artefact is specified per electrode
row_spec_elek = find(ismember(elec_art, (find(~strcmp(dataBase(subj).tb_events.electrodes_involved_onset,'all')))));     % Electrode contains label artefact, though does not contain label all  
art_spec_elek = dataBase(subj).tb_events.electrodes_involved_onset(row_spec_elek);                                       % Find name of electrode involved
art_elek_start = dataBase(subj).tb_events.sample_start(row_spec_elek);
art_elek_stop = dataBase(subj).tb_events.sample_end(row_spec_elek);

% Find sample start stimulation, because artefacts before stimart are
% ignored
stim_start = find(strcmp(dataBase(subj).tb_events.trial_type,'electrical_stimulation'));
start_stim = dataBase(subj).tb_events.sample_start(stim_start(1));

% For PRIOS02 there are many artefacts before stimulation start these are
% ignored.
for i = 1:size(art_spec_elek,1)
    if art_elek_stop(i) > start_stim

        if iscell(art_elek_start)
            art_elek_start = str2double(art_elek_start);
        end

        if iscell(art_elek_stop)
            art_elek_stop = str2double(art_elek_stop);
        end

        ev_artefact_elek = zeros(1,100);        % Preallocation with a guess of the size, structs with empty cells are removed.
        % When there are multiple artefacts, then the times samples are concatenaded
        % in the same array. So ev_artefact is not overwritten
        ev_artefact_elek = [ev_artefact_elek, art_elek_start(i):art_elek_stop(i)];    %#ok<AGROW>

        dataBase(subj).ev_artefact_elek = ev_artefact_elek;

        clear i

    else
        % Do nothing because the artefact is already over before stimulation
        % started
    end
end
%     newStr = extractBetween(art_spec_elek(1,:),',')   ;                 % Split when multiple electrodes are named in one annotaion

   
     
    %% Remove Burst suppression periods 
    BS_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.trial_type,'burst_suppression'));
    BS_stop = dataBase(subj).tb_events.sample_end(strcmp(dataBase(subj).tb_events.trial_type,'burst_suppression'));

    if iscell(BS_start)
        BS_start = str2double(BS_start);
    end
    if iscell(BS_stop)
        BS_stop = str2double(BS_stop);
    end

    BS_sig = zeros(1,100);          % Preallocation with a guess of the size, structs with empty cells are removed.
    for i=1:size(BS_start,1)
        BS_sig = [BS_sig, BS_start(i):BS_stop(i)];  %#ok<AGROW>
    end

    dataBase(subj).Burstsup = BS_sig;

    
     %% Remove seizure periods 
     % This is yet no electrode specific
    SZ_start = dataBase(subj).tb_events.sample_start(strcmp(dataBase(subj).tb_events.trial_type,'seizure'));
    SZ_stop = dataBase(subj).tb_events.sample_end(strcmp(dataBase(subj).tb_events.trial_type,'seizure'));

    if iscell(SZ_start)
        SZ_start = str2double(SZ_start);
    end
    if iscell(SZ_stop)
        SZ_stop = str2double(SZ_stop);
    end

    SZ_sig = zeros(1,100);                  % Preallocation with a guess of the size, structs with empty cells are removed.
    for i=1:size(SZ_start,1)
       SZ_sig  = [SZ_sig, SZ_start(i):SZ_stop(i)];              %#ok<AGROW> % When there are multiple SZ, then the times samples are concatenaded in the same array
    end

    dataBase(subj).Seizure = SZ_sig;

    
    %% Remove stimulus pairs with less than minimal interstim time
    % Stimulations with too little interstimulus time need to be removed
    % since it cannot be guaranteed that the signal went back to baseline
    % before the next stimulus was given.

     idx_tbelec = contains(dataBase.tb_events.trial_type,'electrical_stimulation');
     tb_stim = dataBase.tb_events(idx_tbelec,:);
     
     % Keep stimulations that are at least 3 seconds apart (because of
     % mechanical problems with the stimulator
     idx_keep5s = [true(1); diff([tb_stim.onset])>=3];
     dataBase.tb_events = tb_stim(idx_keep5s,:);
       
         
    %% Unique stimulation pairs
    stimpair = dataBase(subj).tb_events.electrical_stimulation_site(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a')) ;
    
    stimnum = NaN(size(stimpair,1),2);
    for stimp = 1:size(stimpair,1)
        stimchans = strsplit(stimpair{stimp},'-');
        for chan = 1:2
            stimnum(stimp,chan) = find(strcmp(stimchans{chan}, dataBase(subj).ch)==1);
        end
    end
    
    % these are set in config_CCEP
    if strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'yes') % take into account the direction (C1-C2 and C2-C1 separately) and the stimulation current
        stimcur = dataBase(subj).tb_events.electrical_stimulation_current(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a'));
        stimelek = [stimnum stimcur];
        
    elseif strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'no') % take into account only the direction
        stimelek = stimnum;
        
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'yes') % take into account only the stimulation current
        stimcur = dataBase(subj).tb_events.electrical_stimulation_current(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a'));
        stimelek = [sort(stimnum,2) stimcur];
        
    elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'no') % do not take stimulation current or direction into account
        stimelek = sort(stimnum,2);
        
    end
    
    [cc_stimsets_all,~,IC_all] = unique(stimelek,'rows');      
     
    %% Remove stimulus pairs with less than minimum number of stimulations
    n = histcounts(IC_all,'BinMethod','integers');                      
    
    if any(diff(n) ~= 0) % if any pair is stimulated a different amount
        
       if n<minstim
           stimremove = n<minstim;               
           % remove stim pairs in both directions
           remove_elec = cc_stimsets_all(stimremove,:);
           remove_stimp = find(cc_stimsets_all(:,2)==remove_elec(:,2) & cc_stimsets_all(:,1)==remove_elec(:,1) |  cc_stimsets_all(:,2)==remove_elec(:,1) &  cc_stimsets_all(:,1)==remove_elec(:,2));

           warning('%s: stimulation pair(s) are stimulated less than all others, these are removed\n',dataBase(subj).sub_label);

           for i = 1:length(remove_stimp)-1
               stimelek(IC_all==remove_stimp(i),:)= NaN;
               stimelek(IC_all==remove_stimp(i+1),:)= NaN;
           end

           remove(:,1) = isnan(stimelek(:,1));
           stimelek(remove,:) = [];

           [cc_stimsets_all,~,IC_all] = unique(stimelek,'rows');
           n = histcounts(IC_all,'BinMethod','integers');
       end
       warning('%s: a stimulation pair is probably stimulated more than others, often no problem \n',dataBase(subj).sub_label)
    end
    
    
  
    %% Number of stimulations per stimulus pair
    % find the amount of time most stimulus pairs are stimulated and set
    % that as maximal number of stimulus pairs
    max_stim = max(n);                    %max_stim = max(n);               % The max stim per stimulation pair DIRECTION!
    
    if strcmp(cfg.dir_avg,'yes')        % dir_avg = 'yes' analyse the average signal of both the positive and negative stimuli
        
        sort_cc_stimsets = sort(cc_stimsets_all,2);
        [cc_stimsets_avg, ~, IC_avg] = unique(sort_cc_stimsets,'rows');
        
         if 2*length(cc_stimsets_avg) ~= length(cc_stimsets_all)                 % When not all stimulation pairs are stimulated in both directions
             Ncount = find(histcounts(IC_avg,length(cc_stimsets_avg))~=2)';      % stimpairs which are stimulated in one direction            
             
             % Allocation
             remove_rows = zeros(length(Ncount),1);
             for i = 1:length(Ncount)
               remove_rows(i,:) = find(IC_avg==Ncount(i,:));                    % Row in which the 'single' stimpair is located
             end
             
             % Remove rows with 'single' stimpairs
              cc_stimsets_all(remove_rows,:) = [];   
              sort_cc_stimsets(remove_rows,:) = [];
              % recalculate IC_avg
              [cc_stimsets_avg, ~, IC_avg] = unique(sort_cc_stimsets,'rows');  
         end
         
    elseif strcmp(cfg.dir_avg,'no')         % dir_avg = 'no' to analyse all signals separately
        
        cc_stimsets_avg = cc_stimsets_all;
        IC_avg = IC_all;
    end
    
    
    %% Determine stimulation pair names
    % pre-allocation
    cc_stimchans_all = cell(size(cc_stimsets_all,1),2);
    cc_stimchans_avg = cell(size(cc_stimsets_avg,1),2);
    cc_stimpnames_all = cell(1,size(cc_stimsets_all,1));
    cc_stimpnames_avg = cell(1,size(cc_stimsets_avg,1));
    
    % determine stimulus pair names for each stimulus
    for stimp = 1:size(cc_stimsets_all,1)
        for chan = 1:2
            cc_stimchans_all{stimp,chan} = dataBase(subj).ch{cc_stimsets_all(stimp,chan)};
        end
        cc_stimpnames_all{stimp} = [cc_stimchans_all{stimp,1} '-' cc_stimchans_all{stimp,2}];
    end
    
    % determine stimulus pair names for averaged stimulus
    for stimp = 1:size(cc_stimsets_avg,1)
        for chan = 1:2
            cc_stimchans_avg{stimp,chan} = dataBase(subj).ch{cc_stimsets_avg(stimp,chan)};
        end
        cc_stimpnames_avg{stimp} = [cc_stimchans_avg{stimp,1} '-' cc_stimchans_avg{stimp,2}];
    end    
  
    dataBase(subj).cc_stimsets_all = cc_stimsets_all;
    dataBase(subj).cc_stimsets_avg = cc_stimsets_avg;
    dataBase(subj).cc_stimchans_all = cc_stimchans_all;
    dataBase(subj).cc_stimchans_avg = cc_stimchans_avg;
    dataBase(subj).stimpnames_all = cc_stimpnames_all;
    dataBase(subj).stimpnames_avg = cc_stimpnames_avg;
    dataBase(subj).max_stim = max_stim;
    
    %% Select epochs
    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    tt = (1:epoch_length*dataBase(subj).ccep_header.Fs)/dataBase(subj).ccep_header.Fs - epoch_prestim;

    % allocation
    cc_epoch_sorted_all = NaN(size(dataBase(subj).data,1),dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets_all,1),t);
    tt_epoch_sorted_all = NaN(dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets_all,1),t); 

    for elec = 1:size(dataBase(subj).data,1)                    % for all channels 
       for ll = 1:size(dataBase(subj).cc_stimsets_all,1)       % for all epochs with > minimum number of stimuli (minstim)

           eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans_all{ll,1}, '-',dataBase(subj).cc_stimchans_all{ll,2}]));

            if size(eventnum,1) > dataBase(subj).max_stim
                events = dataBase(subj).max_stim;
            else
                events = size(eventnum,1);
            end
                
            for n = 1:events
                if dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1< 0
                     % do nothing, because epoch starts before the start of
                     % a file (and gives an error)
                elseif ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),ev_artefact_all)
                       % do nothing, because part of artefact, therefore
                       % this event is not used
                       
                       % When the concerning electrode is part of a electrode-specified artefact, that is not applicable to the whole electrode-dataset
                elseif ismember(dataBase.ch{elec}, art_spec_elek) && ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),ev_artefact_elek)
                    % When eventnum is part of ev_artefact_elek AND the
                    % specified electrode
                                                                              
                    % do nothing, because part of an electrode specified artefact, therefore
                    % this event is not used

                        
                elseif ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),BS_sig)
                    % do nothing because stimulation is part of burst
                    % suppression period and therefore not thrustworthy
                    
                elseif ismember(dataBase(subj).tb_events.sample_start(eventnum(n)),SZ_sig)
                    % do nothing because stimulation is part of a seizure
                    % seizure periods can therefore not be scored
                    
                else
                    cc_epoch_sorted_all(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                    tt_epoch_sorted_all(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);
                                         
                end
            end
        
       end
       
   end

    %% Average epochs
    
    % preallocation
    cc_epoch_sorted_avg = NaN(size(cc_epoch_sorted_all,1),size(cc_stimsets_avg,1),size(cc_epoch_sorted_all,4)); % [channels x stimuli x samples]
    cc_epoch_sorted_select = NaN(size(cc_epoch_sorted_all,1), size(cc_stimsets_avg,1), max_stim*2, size(cc_epoch_sorted_all,4));   % avg_stim*sum(IC_avg==2)   Wat nergens op slaat trouwens...
    
    for ll = 1:max(IC_avg)                     % Takes every value between 1 and max while some numbers are not used, therefore the next line

            stimps = find(IC_avg == ll);
            
        if any(~isnan(cc_epoch_sorted_all(1,:,stimps(1),1))) && any(~isnan(cc_epoch_sorted_all(1,:,stimps(2),1))) % Check whether the stimpair is stimulated in both directions.      
            
            % Average ALL stimuli given to certain stimpair
            selection = cc_epoch_sorted_all(:,:,IC_avg==ll,:);         
            selection_avg =  squeeze(nanmean(selection,2));
            
            while size(size(selection_avg),2) >2
                selection_avg =  squeeze(nanmean(selection_avg,2));
            end

            cc_epoch_sorted_avg(:,ll,:) = selection_avg;
            cc_epoch_sorted_select(:,ll,:,:) = reshape(selection,size(selection,1), size(selection,2)*size(selection,3) ,size(selection,4));          
                      
        end
    end
    

    %% Find rows which are not stimulated in both directions anymore because
% Either one of the stimulation pairs is removed due to burst suppression
% or artefact
if any(find(isnan(cc_epoch_sorted_avg(1,:,1)))>0)
    
    % Avg
    remove_stimp_avg = find(isnan(cc_epoch_sorted_avg(1,:,1)));
    cc_epoch_sorted_avg(:,remove_stimp_avg,:) = [];

    % Select
    cc_epoch_sorted_select(:,remove_stimp_avg,:,:) = [];

    % All
    remove_stimp_all = zeros(size(remove_stimp_avg,2),2);
    for i = 1:size(remove_stimp_avg,2)
        remove_stimp_all(i,1:2) = find(IC_avg == remove_stimp_avg(i))';  
    end
    remove_stimp_all = unique(reshape(remove_stimp_all,[],1));

    cc_epoch_sorted_all(:,:,remove_stimp_all,:) = [];

    % tt
    tt_epoch_sorted_all(:,remove_stimp_all,:) = [];
    
    
    cc_stimsets_all(remove_stimp_all,:) = [];
    cc_stimsets_avg(remove_stimp_avg,:) = [];
    cc_stimchans_all(remove_stimp_all,:) = [];
    cc_stimchans_avg(remove_stimp_avg,:) = [];
    cc_stimpnames_all(:,remove_stimp_all) = [];
    cc_stimpnames_avg(:,remove_stimp_avg) = [];
    fprintf('...%d are removed because only one direction was stimulated... \n',remove_stimp_avg(1,:))
   
end




    dataBase(subj).cc_stimsets_all = cc_stimsets_all;
    dataBase(subj).cc_stimsets_avg = cc_stimsets_avg;
    dataBase(subj).cc_stimchans_all = cc_stimchans_all;
    dataBase(subj).cc_stimchans_avg = cc_stimchans_avg;
    dataBase(subj).stimpnames_all = cc_stimpnames_all;
    dataBase(subj).stimpnames_avg = cc_stimpnames_avg;
    dataBase(subj).max_stim = max_stim;
    
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted_all;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted_all;
    dataBase(subj).tt = tt;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    dataBase(subj).cc_epoch_sorted_select_avg = cc_epoch_sorted_select;    
    dataBase(subj).stimpnames = cc_stimpnames_all;

    
end
end