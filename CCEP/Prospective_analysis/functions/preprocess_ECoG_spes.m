function dataBase = preprocess_ECoG_spes(dataBase,cfg,avg_stim)

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
    avg_stim = [];
end

for subj = 1:size(dataBase,2)
    
    %% define artefact period - stimuli during an artefact become NaNs
    
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
        dataBase(subj).ev_artefact = [ev_artefact, ev_artefact_start(i):ev_artefact_stop(i)]; 
    end
    
    clear i
    
    %% Remove stimulus pairs with less than minimal interstim time
    % Stimulations with too little interstimulus time need to be removed
    % since it cannot be guaranteed that the signal went back to baseline
    % before the next stimulus was given.

     idx_tbelec = contains(dataBase.tb_events.trial_type,'electrical_stimulation');
     tb_stim = dataBase.tb_events(idx_tbelec,:);
     
     idx_keep5s = [true(1); diff([tb_stim.onset])>=4.9];
     dataBase.tb_events = tb_stim(idx_keep5s,:);
       
     
    %% unique stimulation pairs
    stimpair = dataBase(subj).tb_events.electrical_stimulation_site(contains(dataBase(subj).tb_events.sub_type,'SPES') & ~contains(dataBase(subj).tb_events.electrical_stimulation_site,'n/a')) ;
    
    
    %%% FOR PRIOS04 GAAT HET MIS ZIE TB_EVENTS REGELS 47-52, MOGELIJK IETS
    %%% INBOUWEN WAARBIJ JE AANGEEFT DAT ER MINIMAAL 5 CHARACTERS ZIJN
    stimnum = NaN(size(stimpair,1),2);
    for stimp = 1:size(stimpair,1)
        stimchans = strsplit(stimpair{stimp},'-');
        for chan = 1:2
            stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(subj).ch)==1);
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
     
    %% remove stimulus pairs with less than minimum number of stimulations
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
    max_stim = median(n);               % The max stim per stimulation pair DIRECTION!
    
    if strcmp(cfg.dir_avg,'yes')        % dir_avg = 'yes' analyse the average signal of both the positive and negative stimuli
        
        sort_cc_stimsets = sort(cc_stimsets_all,2);
        [cc_stimsets_avg, ~, IC_avg] = unique(sort_cc_stimsets,'rows');
        
         if 2*length(cc_stimsets_avg) ~= length(cc_stimsets_all)                % When not all stimulation pairs are stimulated in both directions
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
    
    
    %% determine stimulation pair names
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
    
  
   
    %% select epochs
    t = round(epoch_length*dataBase(subj).ccep_header.Fs);
    tt = (1:epoch_length*dataBase(subj).ccep_header.Fs)/dataBase(subj).ccep_header.Fs - epoch_prestim;
    
    %%% HIER MOET IETS INGEBOUWD WORDEN DAT IK KAN KIEZEN WELKE PROPOFOL STIMULATIE
    %%% IK WIL, NU IS HET TELKENS DE LAATSTE WAT NIET ALTIJD CORRECT IS.
    %%% KOMT DOOR MAX STIM = 1, HEB DIT OOK WEL OP 5 GEHAD MAAR DAN KRIJG
    %%% JE VEEL RIJEN MET NaN EN DAT KLOPT OOK NIET.
    
    % allocation
    cc_epoch_sorted_all = NaN(size(dataBase(subj).data,1),dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets_all,1),t);
    tt_epoch_sorted_all = NaN(dataBase(subj).max_stim,size(dataBase(subj).cc_stimsets_all,1),t); % samplenumbers for each epoch % Dit is nu een 2xstimpairXtt matrix
    
    for elec = 1:size(dataBase(subj).data,1)                    % for all channels
        for ll = 1:size(dataBase(subj).cc_stimsets_all,1)       % for all stimulation pair directions with the minimum number of stimuli (minstim)
            if strcmp(cfg.dir,'no')
                
                % Find the stimulationnumbers on which a stimulation pair is stimulated in both directions.
                eventnum1 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans_all{ll,1}, '-',dataBase(subj).cc_stimchans_all{ll,2}])); % Positive stimulation
                eventnum2 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans_all{ll,2}, '-',dataBase(subj).cc_stimchans_all{ll,1}])); % Negative stimulation
                eventnum = [eventnum1;eventnum2];
                
            elseif strcmp(cfg.dir,'yes')
                
                eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans_all{ll,1}, '-',dataBase(subj).cc_stimchans_all{ll,2}]));
            end
            
            if size(eventnum,1) > dataBase(subj).max_stim
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
                    cc_epoch_sorted_all(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
                    tt_epoch_sorted_all(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);
                end
            end
        end
    end
    
    %% average epochs
    if isempty(avg_stim)
        avg_stim = max_stim;
    end
    
    % preallocation
    cc_epoch_sorted_avg = NaN(size(cc_epoch_sorted_all,1),size(cc_stimsets_avg,1),size(cc_epoch_sorted_all,4)); % [channels x stimuli x samples]
    cc_epoch_sorted_select = NaN(size(cc_epoch_sorted_all,1),size(cc_stimsets_avg,1),avg_stim*sum(IC_avg==2),size(cc_epoch_sorted_all,4)); % [channels x stimuli x selected trials x samples[
    
    for ll = 1:max(IC_avg)                     % Takes every value between 1 and max while some numbers are not used, therefore the next line
        if sum(IC_avg==ll)>1                   % Check whether the stimpair is stimulated in both directions.
            
            selection = cc_epoch_sorted_all(:,1:avg_stim,IC_avg==ll,:);
            selection_avg =  squeeze(nanmean(selection,2));

            while size(size(selection_avg),2) >2
                selection_avg =  squeeze(nanmean(selection_avg,2));
            end

            cc_epoch_sorted_avg(:,ll,:) = selection_avg;
            cc_epoch_sorted_select(:,ll,:,:) = reshape(selection,size(selection,1),size(selection,2)*size(selection,3),size(selection,4));
        end
    end
    
    dataBase(subj).cc_epoch_sorted = cc_epoch_sorted_all;
    dataBase(subj).tt_epoch_sorted = tt_epoch_sorted_all;
    dataBase(subj).tt = tt;
    dataBase(subj).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
    dataBase(subj).cc_epoch_sorted_select_avg = cc_epoch_sorted_select;
   
    dataBase(subj).stimpnames = cc_stimpnames_all;
    
    fprintf('...%s has been epoched and averaged... \n',dataBase(subj).sub_label)
    
end
end