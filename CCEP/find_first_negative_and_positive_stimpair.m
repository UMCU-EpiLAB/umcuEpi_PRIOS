%% Determine the tb_events.sample_start voor alleen de eerste van het eventnum
%  cc_epoch_sorted(elec,n,ll,:) = dataBase(subj).data(elec,dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs));
%  tt_epoch_sorted(n,ll,:) = dataBase(subj).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).ccep_header.Fs)+1:dataBase(subj).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).ccep_header.Fs);

% Use eventunum and cfg.dir as described below
% if strcmp(cfg.dir,'no')
%     eventnum1 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}])); % Positive stimulation
%     eventnum2 = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,2}, '-',dataBase(subj).cc_stimchans{ll,1}])); % Negative stimulation
%     eventnum = [eventnum1;eventnum2];
% elseif strcmp(cfg.dir,'yes')
%     eventnum = find(strcmp(dataBase(subj).tb_events.electrical_stimulation_site,[dataBase(subj).cc_stimchans{ll,1}, '-',dataBase(subj).cc_stimchans{ll,2}]));
% end          

