function [dataBase] = detect_n1peak_ECoG_ccep(dataBase, cfg)

% Function for detecting the the N1 peaks of CCEPs

% INPUTS:
% - database
%   Structure with metadata and averaged epoched data split into
%   [channels x stimulations x time].
% - amplitude_thresh
%   Threshold factor that needs to be exceeded to detect a peak as
%   significant peak. This factor is multiplies the minimal pre simulation
%   standard deviation of 50 uV. amplitude_thresh of 3.4 recommended for
%   conservative algorithm (50 uV * 3.4 = 170 uV).
% - n1_peak_range
%   End point (in ms) of the range in which the algorithm will search for the
%   peak of the N1. Algorithm will search between 10ms and the end point.
%   End point of 90 ms recommended for conservative algorithm (80ms has
%   similar performance.

% OUTPUTS:
% - database
%   To the database structure the following matrices will be added:
% - n1_peak_sample
%   matrix[channels x stimulations] containing latency of the significant
%   n1 peaks that are detected.
% - n1_peak_amplitude
%   matrix[channels x stimulations] containing the amplitude of the
%   significant n1 peaks that are detected.

% Peakfinder function is used for peak detection, add to path

% This function is heavily based on the function and algorithm developed by
% Dorien van Blooijs during her master's thesis 'Improving the SPES protocol
% by automating ER and DR detection and evaluation of the spatial relation
% between ERs and DRs. Algorithm is validated during the thesis and used in
% publication (van Blooijs et al., 2018).

% Some slight differences:
% - While algorithm of van Blooijs detects both positive as negative peaks,
%   this algorithm only selects negative peaks
% - After validation the algorithm of van Blooijs used 125uV as amplitude threshold
%   and 100ms af n1 peak range, this algorithm performes best with 140 uV
%   and 80 ms. A reason for this possibly is the other ages of patients which are
%   used during validation. The first validation is done is only younger
%   (9yrs - 12 yrs.) patients, while later validation is also done in older
%   patients. Age seems to effect the characteristics of CCEPs and
%   therefore also the detection.
% - look at peak range....

% This version of the algorithm is validated in three patients (age 9, 21 and 50)
% Validation done by comparing the performance of the code and visual
% assesment. Algorithm optimized by setting different parameters and
% comparing their performances. (see parameters optimize function).

% For a conservative algorithm (high specificity of at least 95%), the following
% parameters are advised: (see validation matrices for performances with
% other parameters e.g. if you want a very sensitive algorithm)
%<<<<<<< HEAD
% - amplitude threshold of 140 uV (minSD * threshold = 50 uV * 2.8)
%   recommended
% - N1 peak range of (10 to) 70 ms is recommended
%=======
% - amplitude threshold of 170 uV (minSD * threshold = 50 uV * 3.4)
%   recommended
% - N1 peak range of (10 to) 90 ms is recommended
%>>>>>>> upstream/master

% FIXED PARAMETERS (that are validated by both van Blooijs & van der Aar):
% - sel = 20, which is how many samples around peak not considered as another peak
% - minSD = 50, minimum standard deviation (in uV) of prestimulus baseline

% original author: Dorien van Blooijs, UMC Utrecht, January 2018
% modified by: Jaap van der Aar, Dora Hermes, Dorien van Blooijs, Giulio Castegnaro, UMC Utrecht, 2019


amplitude_thresh = cfg.amplitude_thresh;
n1_peak_range = cfg.n1_peak_range;
epoch_prestim = cfg.epoch_prestim;
epoch_length = cfg.epoch_length;

bad_channels = find(contains(dataBase.tb_channels.status,'bad')==1);

%% Script
% iterate over all subjects in database
for subj = 1:length(dataBase)
    % iterate over all their runs
    %     for runs = 1:length(dataBase(subj).metadata)
    
    % output in channels X stimulations X [latency amplitude]
    n1_peak = NaN(size(dataBase(subj).cc_epoch_sorted_avg,1), ...
        size(dataBase(subj).cc_epoch_sorted_avg,2),2);
    
    % for every averaged stimulation
    for jj = 1:size(dataBase(subj).stimpnames_avg,2)                         % for every averaged stimulation
        
        % for every channel
        for ii = 1:size(dataBase(subj).ch,1)                                 % for every channel
            
            % create time struct
            tt = (1:epoch_length*dataBase(subj).ccep_header.Fs) / ...
                dataBase(subj).ccep_header.Fs - epoch_prestim;
            
            % baseline subtraction: take median of part of the averaged signal for
            % this stimulation pair before stimulation, which is the half of the
            % epoch
            baseline_tt = tt>-2 & tt<-.3;                   %%%% THIS WAS UNTIL 0.1%%%%
            
            signal_median = median(dataBase(subj).cc_epoch_sorted_avg(ii,jj,baseline_tt),3);
            
            % subtract median baseline from signal
            new_signal = squeeze(dataBase(subj).cc_epoch_sorted_avg(ii,jj,:)) - signal_median;
            % testplot new signal: plot(tt,squeeze(new_signal))
            
            % take area before the stimulation of the new signal and calculate its SD
            pre_stim_sd = std(new_signal(baseline_tt));
            
            % if the pre_stim_sd is smaller that the minimally needed SD,
            % which is validated as 50 uV, use this the minSD as pre_stim_sd
            if pre_stim_sd < 50
                pre_stim_sd = 50;
            end
            
            
           if ismember(dataBase.task_name,'task-SPESprop')
                
               Fs = dataBase.ccep_header.Fs;
               Fn = Fs/2;
               
             % Take a small part of the signal before the stimulation artefact to determine the frequencies
%              baseline_tt = tt>-1.0 & tt<-.1;
%              baseline_tt_plot = tt(baseline_tt ==1);
%              new_signal_part = new_signal(baseline_tt);
%              norm_part = sqrt(new_signal_part.^2);
%              [b,a] = butter(2,[32/Fn 39/Fn],'stop');
%              
%              norm_part_filt = filtfilt(b,a,norm_part);
%              [pwwP,fP] = periodogram(norm_part, rectwin(length(norm_part)), [], Fs);                                % without filter           
%              [pwwP30Hz,fP30Hz] = periodogram(norm_part_filt, rectwin(length(norm_part_filt)), [], Fs);              % with filter           
% 
%              figure()
%              subplot(2,1,1)
%              plot(fP,pwwP);
%              ylabel('PSD'); xlabel('Frequency (Hz)');
%              xlim([0 80]);    
%              ylim([0 200])
%              title('Part of signal before stim artefact')
%              
%              subplot(2,1,2)
%              plot(fP30Hz,pwwP30Hz);
%              ylabel('PSD'); xlabel('Frequency (Hz)');
%              xlim([0 80]);    
%              ylim([0 200])
%              title('Part of signal before stim artefact, 30 Hz filtered')
%               
%             % check
%             figure()
%             subplot(2,1,1)
%             plot(tt, new_signal_stimart)
%             title('Original, without stimulation artefact')
%             ylim([-600 600]); xlim([-1 -0.1])
%             
%             subplot(2,1,2)
%             plot(baseline_tt_plot, norm_part_filt)
%             title('Filtered 36 Hz')
%             ylim([-600 600])
%             
%               
              
            % Remove period of the stimulation artefact
            % replace with zero 
            new_signal_stimart = new_signal;
            new_signal_stimart(tt>-0.005 & tt<0.005) = 0;


            % 50 Hz filter combined with a bandstop filter between 33 and
            % 38 Hz, and a bandpass filter for all EEG-like frequencies
%             signal_norm = sqrt(new_signal_stimart.^2);
            
            % Butterworth
%              [b50,a50] = butter(2,[49/Fn 50/Fn],'stop');
%              [b36,a36] = butter(2,[33/Fn 38/Fn],'stop');
% %              
             [bP,aP] = butter(2,[0.1, 33]/Fn,'bandpass');
           
            % concatenate two filters for two known incorrect frequencies
%             new_signal_filt = filtfilt(b50,a50,signal_norm);
%             new_signal_filt_B = filtfilt(b36,a36,new_signal_filt);
%             
            % band pass filter for EEG valid frequencies
            new_signal_filt_pass = filtfilt(bP,aP, new_signal_stimart);
            
% %             new_signal__withoutNorm = filtfilt(b,a,new_signal_stimart);
%             
% %             [pww,f] = periodogram(new_signal,rectwin(length(new_signal)),[],Fs);                                % with stimulation artefact, without filter           
% %             [pwwStimA,fStimA] = periodogram(new_signal_stimart,rectwin(length(new_signal_stimart)),[],Fs);    % without stimulation artefact, without filter
% %              [pwwAF,fAF] = periodogram(new_signal_filt_ori,rectwin(length(new_signal_filt_ori)),[],Fs);        % with stimulation artefact, with filter 
%             [pwwF,fF] = periodogram(new_signal_filt_B,rectwin(length(new_signal_filt_B)),[],Fs);                    % without stimulation artefact, with filter
%             [pwwP,fP] = periodogram(new_signal_filt_pass, [], [], Fs);
%             
% %             [pwwWN,fWN] = periodogram(new_signal__withoutNorm,rectwin(length(new_signal__withoutNorm)),[],Fs);                    % without stimulation artefact, with filter

%             
%             figure() 
%             subplot(2,1,1)
%             plot(fF,pwwF);
%             ylabel('PSD'); xlabel('Frequency (Hz)');
%             xlim([2 60]); ylim([0 200])
%             title('Without stimulation artefact and with 50 and 36 Hz filter')
%             
%             subplot(2,1,2)
%             plot(fP,pwwP);
%             ylabel('PSD'); xlabel('Frequency (Hz)');
%             xlim([2 60]); ylim([0 200])
%             title('Without stimulation artefact and with bandpass filter')
% 
%             
%              % Check
%             figure()
%             subplot(3,1,1)
%             plot(tt, new_signal_stimart)
%             title('without stimulation artefact, not filtered')
%             ylim([-600 600]); xlim([-0.2 1.5])
%             
%             subplot(3,1,2)
%             plot(tt,new_signal_filt_B)
%             title('50 & 36 Hz filtered signal')
%             ylim([-600 600]); xlim([-0.2 1.5])
%          
%             subplot(3,1,3)
%           % Check
%             figure(1)
%             subplot(3,1,1)
%             plot(tt,new_signal)
%             title(sprintf('original %s, %s',dataBase.stimpnames_avg{jj},dataBase.ch{ii}))
% %             ylim([-600 600])
%             subplot(3,1,2)
%             plot(tt, new_signal_stimart)
%             title('without stimulation artefact, zeros')
% %             ylim([-600 600])
%             
%            subplot(3,1,3)           
%             plot(tt,new_signal_filt_pass)
%             title('bandpass filtered signal')
% %             ylim([-300 1000]); 
%             xlim([-1.2 1.4])
%             
%             figure(3)
%             plot(tt,new_signal_stimart,tt,new_signal_filt_pass,'LineWidth',5)
%             legend('original','filtered')
% %          
              new_signal = new_signal_filt_pass ;
           end
        
           
     
            % Place NaN when the electrode is stimulated
            if  ii == dataBase(subj).cc_stimsets_avg(jj,1) || ...
                    ii == dataBase(subj).cc_stimsets_avg(jj,2)
                n1_peak_sample = NaN;
                n1_peak_amplitude = NaN;
                
            elseif ismember(ii,bad_channels) % place NaN when electrode is a bad channel
                n1_peak_sample = NaN;
                n1_peak_amplitude = NaN;

            else
                
                % use peakfinder to find all positive and negative peaks and their
                % amplitude.
                % tt are the samples of the epoch based on the Fs and -2.5 to 2.5
                % seconds sample of the total epoch
                % As tt use first sample after timepoint 0
                % till first sample after 0,5 seconds (rougly 1000 samples)
                % sel = 20 , which is how many samples around a peak not considered as another peak
                % ccep_peakfinder(x0, sel, thresh, extrema)
                [all_sampneg, all_amplneg] = ccep_peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),20,[],-1); % yet no threshold
                
                % If the first selected sample is a peak, this is not a real peak,
                % so delete
                all_amplneg(all_sampneg==1) = [];
                all_sampneg(all_sampneg==1) = [];
                
                % convert back timepoints based on tt, substract 1 because before
                % the first sample after stimulation is taken
                all_sampneg = all_sampneg + find(tt>0,1) - 1;
                
                % set the starting range in which the algorithm looks
                % for a peak. At least 18 samples are necessary because
                % of the micromed amplifier does not record the
                % stimulated electrode before this. Peak detection
                % start 9 ms after stimulation, which is 19 samples
                % after stimulation
                n1_samples_start = find(tt>0.009,1);
                
                % find first sample that corresponds with the given n1
                % peak range
                n1_samples_end = find(tt>(n1_peak_range/1000),1);
                
                
                % for N1, first select the range in which the N1 could appear, and
                % select the peaks found in this range
                temp_n1_peaks_samp = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
                temp_n1_peaks_ampl = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
                
                % if peak(s) found, select biggest peak
                if ~isempty(temp_n1_peaks_samp)
                    max_n1_ampl = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
                    n1_peak_sample = temp_n1_peaks_samp(max_n1_ampl(1));
                    n1_peak_amplitude = temp_n1_peaks_ampl(max_n1_ampl(1));
                    % otherwise give the amplitude the value NaN
                elseif isempty(temp_n1_peaks_samp)
                    n1_peak_amplitude = NaN;
                end
                
                % if N1 exceeds positive threshold, it is deleted
                if temp_n1_peaks_ampl > 0
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
                
                % when peak amplitude is saturated, it is deleted
                if abs(n1_peak_amplitude) > 3000
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
                
                % if the peak is not big enough to consider as a peak, assign NaN
                if abs(n1_peak_amplitude) < amplitude_thresh* abs(pre_stim_sd)
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                end
            end
            
            
            % add properties to output frame
            n1_peak(ii,jj,1) = n1_peak_sample;
            n1_peak(ii,jj,2) = n1_peak_amplitude;
            
        end
    end
    
    % write n1_peak (sample and amplitude) to database
    dataBase(subj).ccep.n1_peak_sample = n1_peak(:,:,1);
    dataBase(subj).ccep.n1_peak_amplitude = n1_peak(:,:,2);
    dataBase(subj).ccep.amplitude_thresh = amplitude_thresh;
    dataBase(subj).ccep.n1_peak_range = n1_peak_range;

end
end