function dataBase = filter_bedArt(dataBase)
% Filter the data to remove noise
% Filtering is preceded by removal of the stimulation artefact. This is
% done with interpolation with the points around the stimulation artefact.
% The stimulation artefact interval is [-1.5 ms : 9 ms] interval which is
% caused by the saturation effect induced by the amplifier.

Fs = dataBase(1).ccep_header.Fs;
Fn = Fs/2;

% Filters: Butterworth, 4th order
[b50,a50] =butter(4,[49/Fn 51/Fn],'stop');
[b36,a36] = butter(4,[33/Fn 38/Fn],'stop');
[b110,a110] = butter(4,[108/Fn 112/Fn], 'stop');
[bP,aP] = butter(4,120/Fn,'low');
 
% fvtool(b36,a36)               % Check the effect of the filter
% fvtool(bP,aP)

% Remove all stimulation artefacts. 
% Interpolate for the stimulation artefact  
for i = 1:size(dataBase,2)
    data = dataBase(i).raw_data;
    
    for event = 1:size(dataBase(i).tb_events,1)
         if strcmp(dataBase(i).tb_events.trial_type(event), 'electrical_stimulation')
             
            % The [-1.5 ms : 9 ms] interval surrounding the stimulus was removed using interpolation
            % to avoid incorrect filtering becuase of the saturation effect induced by the amplifier 
            % (the stimulation artefact) @VanBlooijs2015, @VantKlooster2011
            
            stimart_start = dataBase(i).tb_events.sample_start(event) - 3;       
            stimart_stop =  dataBase(i).tb_events.sample_start(event) +19;       
            
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
         end        
    end
          
%   % plot the whole signal without stimulation artefacts
%     figure()
%     ch = 5;
%     subplot(2,1,1)
%     plot(dataBase(i).raw_data(ch,:))
%     legend('original')
%     subplot(2,1,2)
%     plot(data(ch,:))
%     legend('without artefact')
%     title('Whole signal with and without the stimulation artefacts')
%     xlabel('time (samples')
     
     
    
    % Filter every signal
    % Preallocation
    signal = zeros(size(data,1), size(data,2));
    signal_filt_36 = zeros(size(data,1), size(data,2));
    signal_filt_3650 = zeros(size(data,1), size(data,2));
    signal_filt_3650110 = zeros(size(data,1), size(data,2));
    signal_filt_pass = zeros(size(data,1), size(data,2));
    
    for ch = 1:size(data,1)   
        signal(ch,:) = data(ch,:);                                          % use the data without stimulation artefact.
        signal_filt_36(ch,:) = filtfilt(b36,a36,signal(ch,:));              % First filter the 36 Hz artefact out
        signal_filt_3650(ch,:) = filtfilt(b50,a50,signal_filt_36(ch,:));
        signal_filt_3650110(ch,:) = filtfilt(b110,a110,signal_filt_3650(ch,:));
        signal_filt_pass(ch,:) = filtfilt(bP,aP, signal_filt_3650110(ch,:));   % Then use the 120 Hz low pass filter to remove all frequencies above 120 Hz.
              
    end

    
% % Periodograms to determine the frequencies in the whole signal
%     t_start = 1704324;
%     t_stop = 2553463;
% 
%     [pww,f] = periodogram(signal(ch,((t_start) : (t_stop))),[],[],Fs);                     % original with bed art
%     [pww_filt,f_filt] = periodogram(signal_filt_pass(ch,((t_start) : (t_stop))), [], [], Fs);      % Filtered signal    
% 
%     % Plot the frequencies
%     figure('Position',[318,133,1222,904])
%     subplot(2,2,1)
%     plot(f,pww,'LineWidth',2)
%     legend('Original')
%     xlim([12 122])
%     xlabel('Frequencies (Hz)')
%     title(sprintf('%s, %s, channel: %s',dataBase(i).sub_label, dataBase(i).run_label,dataBase(i).ch{ch}))
% 
%     subplot(2,2,3)
%     plot(f_filt ,pww_filt,'LineWidth',2)
%     legend('Filtered')
% %     ylim([0 300])
%     xlim([12 122])
%     title(sprintf('%s, %s, channel: %s',dataBase(i).sub_label, dataBase(i).run_label,dataBase(i).ch{ch}))
%     xlabel('Frequencies (Hz)')
% 
%     % Plot the signal response filtered and original  
%     subplot(2,2,[2 4])
%     plot(signal_filt_pass(ch,(t_start : t_start+500)));
%     hold on
%     plot(signal(ch,(t_start : t_start+500)));
% %     ylim([-1000 2000])
%     legend('Filtered','Original')
%     title(sprintf('%s, %s, channel: %s',dataBase(i).sub_label, dataBase(i).run_label,dataBase(i).ch{ch}))
%     
%     % Save figure
%     Frequencies = fullfile(myDataPath.CCEPpath,dataBase(i).sub_label,'/',...
%     [dataBase(i).sub_label '_' dataBase(i).run_label '_frequencies_withoutBedArt' ])
%     set(gcf,'PaperPositionMode','auto');
%     print('-dpng','-r300',Frequencies);
%     savefig(Frequencies);

      
dataBase(i).data = signal_filt_pass;    
end


    

% subtract the median .
end