function dataBase = filter_bedArt(dataBase)
% Filter the data to remove noise
% Filtering is preceded by removal of the stimulation artefact. This is
% done with interpolation with the points around the stimulation artefact.
% The stimulation artefact interval is [-1.5 ms : 9 ms] interval which is
% caused by the saturation effect induced by the amplifier.

Fs = dataBase(1).ccep_header.Fs;
Fn = Fs/2;

% Filters: Butterworth, 4th order
[b50,a50] = butter(4,[49/Fn 51/Fn],'stop'); % band stop 49-51Hz
[b36,a36] = butter(4,[33/Fn 38/Fn],'stop'); % band stop 33-38Hz
[b110,a110] = butter(4,[108/Fn 112/Fn], 'stop'); % band stop 108-112 Hz
[bP,aP] = butter(4,120/Fn,'low'); % lowpass filter untill 120 Hz
 
interplwindow = 40;
% fvtool(b36,a36)               % Check the effect of the filter
% fvtool(bP,aP)

% Remove all stimulation artefacts. 
% Interpolate for the stimulation artefact  
for i = 1:size(dataBase,2)
    data = dataBase(i).raw_data;
    
    for event = 1:size(dataBase(i).tb_events,1)
         if strcmp(dataBase(i).tb_events.trial_type(event), 'electrical_stimulation')
             
            % The [-1.5 ms : 9 ms] = [-3 samples: 19 samples] interval surrounding the stimulus was removed using interpolation
            % to avoid incorrect filtering becuase of the saturation effect induced by the amplifier 
            % (the stimulation artefact) @VanBlooijs2015, @VantKlooster2011
            
            stimart_start = dataBase(i).tb_events.sample_start(event) - 3;       
            stimart_stop =  dataBase(i).tb_events.sample_start(event) +19;       
            
            % Remove the stimulation artefact in every channel
            for channel = 1:size(data,1)                                   
                 data(channel,(stimart_start:stimart_stop)) = NaN;            % Make the period of the stimulation artefact as NaN.
                
                 % interpolate between the points of the stimulation artefact.
                 % Use 40 samples before and after the artefact to calculate the value in the [-1.5 ms : 9 ms] interval. 
                 Interpol_period = data(channel,((stimart_start - interplwindow) : (stimart_stop + interplwindow))) ;              % find the signal 40 samples before and after stimulation artefact,
                 IX = 1:numel(Interpol_period);
                 tf = isnan(Interpol_period);                                                        % Find the NaN's (interval that has to be interpolated)
                 Interpol_period(tf) = interp1(IX(~tf),Interpol_period(~tf),IX(tf));                 % Interpolate between the value -40 and +40 samples around the stimulation artefact
%                  Interpol_period2(tf) = interp1(IX(~tf),Interpol_period(~tf),IX(tf),'spline');                 
                 % spline might be more interesting for future analysis
                 
                 data(channel, ((stimart_start - interplwindow) : (stimart_stop + interplwindow))) = Interpol_period;
                       
                
            end
            
         end        
    end
          
    
               
                 
% Filter every signal
% Preallocation
signal = data'; % [samples x channels]
signal_filt_36 = filtfilt(b36,a36,signal);              % First filter the 36 Hz artefact out
signal_filt_3650 = filtfilt(b50,a50,signal_filt_36);
signal_filt_3650110 = filtfilt(b110,a110,signal_filt_3650);
signal_filt_pass = filtfilt(bP,aP, signal_filt_3650110);   % Then use the 120 Hz low pass filter to remove all frequencies above 120 Hz.
    
      
dataBase(i).data = signal_filt_pass';    
end 

end