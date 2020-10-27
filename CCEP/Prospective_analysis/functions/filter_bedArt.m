function data_filt = filter_bedArt(dataBase, cfg)

% Filter
Fs = dataBase(1).ccep_header.Fs;
Fn = Fs/2;
epoch_length = cfg.epoch_length;         
epoch_prestim = cfg.epoch_prestim;
tt = (1:epoch_length*Fs) / Fs - epoch_prestim;



% Remove all stimulation artefacts. 
% Interpolate for the stimulation artefact  
for i = 1:size(dataBase,2)
    data = dataBase(i).raw_data;
    
    for event = 1:size(dataBase(i).tb_events,1)
         if strcmp(dataBase(i).tb_events.trial_type(event), 'electrical_stimulation')
            stimart_start = dataBase(i).tb_events.sample_start(event) - 3;
            stimart_stop =  dataBase(i).tb_events.sample_end(event) +12;
 
            % Remove the stimulation artefact in every channel
            for channel = 1:size(data,1)                                   
                 data(channel,(stimart_start:stimart_stop)) = NaN;            % Make the period of the stimulation artefact as NaN.
                
                 % interpolate between the points of the stimulation artefact.
                 % Use 40 samples before and after the artefact. 
                 Interpol_period = data(channel,((stimart_start - 40) : (stimart_stop + 40))) ;              % find the signal 40 samples before and after stimulation artefact,
                 IX = 1:numel(Interpol_period);
                 tf = isnan(Interpol_period);                                                        % Find the NaN's 
                 Interpol_period(tf) = interp1(IX(~tf),Interpol_period(~tf),IX(tf));                 % Interpolate between the -40 and +40 samples around the stimulation artefact
             
                 data(channel, ((stimart_start - 40) : (stimart_stop + 40))) = Interpol_period;
                        
            end
            
             % Plot the event in all channels
            for Plotje = 1:size(data,1)
                figure('Position',[852,636,560,420])
                plot(dataBase(i).raw_data(Plotje,((stimart_start - 40) : (stimart_stop + 200))),'LineWidth',2) 

                hold on
                plot(data(Plotje,(stimart_start-40:stimart_stop+200)),'LineWidth',2)
                title(sprintf('SPESclin, PRIOS03, %s, %s',dataBase(i).ch{Plotje}, dataBase(i).tb_events.electrical_stimulation_site{event}))
                legend('Original','Interpolation','Location','southeast')
                xlabel('Samples')
                ylabel('Event per electrodes (µV)')

            end



         end
         
    end
    
    
end

    

% Filters: Butterworth, 4th order
[b36,a36] = butter(4,[33/Fn 38/Fn],'stop');
[bP,aP] = butter(4,120/Fn,'low');
% concatenate two filters for two known incorrect frequencies
new_signal_filt_36= filtfilt(b36,a36,new_signal_stimart);

% band pass filter for EEG valid frequencies
new_signal_filt_pass = filtfilt(bP,aP, new_signal_filt_36);

figure()
plot(tt,new_signal_stimart,tt,new_signal_filt_pass,'LineWidth',1)
legend('original','filtered')
title('stop: 36 Hz, pass: <120Hz')

% subtract the median .
end