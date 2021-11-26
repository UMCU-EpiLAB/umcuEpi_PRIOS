function P1_latency(dataBase_clin, dataBase_prop,cfg, myDataPath)
% Determine the amplitude and latency of the P1 and the highest point before N1
% Necessary to determine the rise and fall times of the N1.

% Script returns with a file with the saved P1 latencies (in
% samples)

% Sampling frequency to transpose the samples to seconds
fs = dataBase_clin.ccep_header.Fs;

% Load the CCEP_checked files
filename = dir([myDataPath.CCEP_allpat,dataBase_clin.sub_label,'*_CCEP_clin_reref_check.mat']);
if ~isempty(filename)
   if exist(fullfile(filename.folder,filename.name),'file')
       ER_checked_clin = load(fullfile(filename.folder,filename.name));
   else
       warning('No N1-checked file is found')
   end
else
    warning('No N1-checked file is found')
end

ER_clin = ER_checked_clin.ccep_clin.n1_peak_sample_check;       
         
% Load the CCEP_checked files
filename = dir([myDataPath.CCEP_allpat,dataBase_prop.sub_label,'*_CCEP_prop_reref_check.mat']);
if ~isempty(filename)
    if exist(fullfile(filename.folder,filename.name),'file')
        ER_checked_prop = load(fullfile(filename.folder,filename.name));
    else
        warning('No N1-checked file is found')
    end
else
    warning('No N1-checked file is found')
end

ER_prop = ER_checked_prop.ccep_prop.n1_peak_sample_check;
 

mode = {'SPES_clin','SPES_prop'};
       
for i = 1:size(mode,2)
    
     if strcmp(mode{i},'SPES_clin')
        dataBase = dataBase_clin;
        
    elseif strcmp(mode{i},'SPES_prop')
        dataBase = dataBase_prop;
        
     end
        
    % Clin and Prop must both have an N1
    r = 0;
    latency_P1.(mode{i}) = NaN(numel(find(~isnan(ER_clin(~isnan(ER_prop))))),1);

    for stim = 1:size(dataBase.stimpnames_avg,2)                          % For each stimpair
        for chan = 1:size(dataBase.ch,1)                                   % For each electrode

        % When both clinical SPES and propofol SPES show an ER
          if ~isnan(ER_clin(chan, stim)) &&  ~isnan(ER_prop(chan, stim))           

              if strcmp(mode{i},'SPES_clin')
                  xN1 = ER_checked_clin.ccep_clin.n1_peak_sample_check(chan,stim)-(cfg.epoch_prestim*fs) ;      % Minus the period before stimulation, in SAMPLES
              elseif strcmp(mode{i},'SPES_prop')
                  xN1 = ER_checked_prop.ccep_prop.n1_peak_sample_check(chan,stim)-(cfg.epoch_prestim*fs) ;      % Minus the period before stimulation
              end

             r = r+1;
             % The N1 peak must at least be 3 samples further than the 9 ms and before 500 ms 
             
             max_latP1 = 1024;     %before 500 ms.
             
             if xN1 > 19+3  && xN1+cfg.epoch_prestim*fs < (max_latP1+cfg.epoch_prestim*fs)-3      
                    
                
                %%%%%%%% P1 %%%%%%%%%%%%%%%%%%%%
                % Find highest peak after the N1 (P1)
                % The interval in wich de P1 is found is [N1 location:500 ms (2*fs+1024 samples = 5120)]

                %%%% should be nice to make this a peak with at least an amplitude of X
%%%% times bigger than the N1 amplitude.

                [~,xP1] = findpeaks(squeeze(dataBase.cc_epoch_sorted_avg(chan,stim,:)));

                % When multiple peaks are found, select the first peak (x)

                if size(xP1,1) > 1    

                   remove = xP1<(xN1+cfg.epoch_prestim*fs) ;                    % Remove all peaks before xN1
                   xP1(remove) = [];
                   remove2 = xP1>((xN1+cfg.epoch_prestim*fs)+max_latP1);             % remove all peaks after 500 ms
                   xP1(remove2) = [];
                   
                   [xP1,~] = (min(xP1) );                       % Select the first peak to be the P1
                end

                % When no peak is found, take the point at 100 ms
                if isempty(xP1)
                    warning('geen xP1')

                end
                
                % Convert to samples after the stimulation artefact
                latency_P1.(mode{i})(r,:) = xP1-(cfg.epoch_prestim*fs);

            end
          
          end
                    
        end
       
    end
    
end


% Find where latency P1 is NaN for SPES-clin and SPES-prop, these are
% skipped in the for loop.
if sum(isnan(latency_P1.SPES_clin)) > 0
    rem_row = find(isnan(latency_P1.SPES_clin));
   
   latency_P1.SPES_clin(rem_row) = [];     
   latency_P1.SPES_prop(rem_row) = [];
  
elseif  sum(isnan(latency_P1.SPES_prop)) > 0
    rem_row = find(isnan(latency_P1.SPES_prop));
   
   latency_P1.SPES_clin(rem_row) = [];   
   latency_P1.SPES_prop(rem_row) = [];
 
end

latency_P1.percentiles_P1_clin = prctile(latency_P1.SPES_clin,[25 50 75]);
latency_P1.percentiles_P1_prop = prctile(latency_P1.SPES_prop,[25 50 75]);

% Save N1 rise and fall times for mathematical analysis
targetFolder = fullfile(myDataPath.CCEPpath, 'Visualise_agreement/N1_compare/P1_latency/');
fileName = [dataBase(1).sub_label,'_P1_latency.mat'];    

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

save([targetFolder,fileName], 'latency_P1');


end


