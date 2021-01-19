function N1_rise_fall(dataBase_clin, dataBase_prop,cfg, myDataPath)
% Determine the amplitude and latency of the P1 and the highest point before N1
% Necessary to determine the rise and fall times of the N1.

% Sampling frequency to transpose the samples to seconds
fs = dataBase_clin.ccep_header.Fs;

% Load the CCEP_checked files
filename = dir([myDataPath.CCEP_allpat,dataBase_clin.sub_label,'*_CCEP_clin_filt_check.mat']);
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
filename = dir([myDataPath.CCEP_allpat,dataBase_prop.sub_label,'*_CCEP_prop_filt_check.mat']);
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
    slope_fall.(mode{i}) = zeros(numel(find(~isnan(ER_clin(~isnan(ER_prop))))),3);                 % preallocation to the size of ER in both clin and prop
    slope_rise.(mode{i}) = zeros(numel(find(~isnan(ER_clin(~isnan(ER_prop))))),3);                 % preallocation to the size of ER in both clin and prop

    for stim = 1:size(dataBase.stimpnames_avg,2)                          % For each stimpair
        for chan = 1:size(dataBase.ch,1)                                   % For each electrode

        % When both clinical SPES and propofol SPES show an ER
          if ~isnan(ER_clin(chan, stim)) &&  ~isnan(ER_prop(chan, stim))           

              if strcmp(mode{i},'SPES_clin')
                  yN1 = ER_checked_clin.ccep_clin.n1_peak_amplitude_check(chan,stim);
                  xN1 = ER_checked_clin.ccep_clin.n1_peak_sample_check(chan,stim)-(cfg.epoch_prestim*fs) ; 
              elseif strcmp(mode{i},'SPES_prop')
                  yN1 = ER_checked_prop.ccep_prop.n1_peak_amplitude_check(chan,stim);
                  xN1 = ER_checked_prop.ccep_prop.n1_peak_sample_check(chan,stim)-(cfg.epoch_prestim*fs) ; 
              end

             r = r+1;
             % The N1 peak must at least be 3 samples further than the 9 ms and before 100 ms 
             if xN1 > 19+3  && xN1+cfg.epoch_prestim*fs <4301-3      
                    
                %%%%%%%% P0 %%%%%%%%%%%%%%%%%%%%          
                % Find highest peak before N1 (P0)
                % The interval in wich de P0 is found is [9 ms (4115 samples) : % N1 location]
                % 4115 = 2*fs(2048) + 19 samples (stim artefact)                    
                [yP0,xP0] = findpeaks(squeeze(dataBase.cc_epoch_sorted_avg(chan,stim,(4115: (xN1+cfg.epoch_prestim*fs)))));

                % When multiple peaks are found, select the highest (y)
                if size(yP0,1) > 1    
                   [yP0,y0_loc] = (max(yP0) );
                   xP0 = xP0(y0_loc);
                end

                % When no peak is found, take the point at 9 ms
                if isempty(xP0)
                    yP0 = dataBase.cc_epoch_sorted_avg(chan,stim,cfg.epoch_prestim*fs+19);
                    xP0 = 0;
                end

                %%%%%%%% P1 %%%%%%%%%%%%%%%%%%%%
                 % Find highest peak before N1 (P1)
                % The interval in wich de P1 is found is [N1 location:100 ms (2*fs+205 samples = 4301)]
                % 4115 = 2*fs(2048) + 19 samples (stim artefact)
                [yP1,xP1] = findpeaks(squeeze(dataBase.cc_epoch_sorted_avg(chan,stim,((xN1+2*fs):4301))));

                % When multiple peaks are found, select the first peak (x)
                % or highest (y)???
                if size(xP1,1) > 1    
                   [xP1,x1_loc] = (min(xP1) );
                   yP1 = yP1(x1_loc);
                end
                
                % When no peak is found, take the point at 100 ms
                if isempty(xP1)
                    yP1 = dataBase.cc_epoch_sorted_avg(chan,stim,2*fs+205);
                    xP1 = 100;
                end

                %%% Determine the fall time (P0 to N1)
                deltaY_fall = yP0-yN1;
                deltaX_fall = xP0 - xN1;
                slope_fall.(mode{i})(r,1) = deltaY_fall / deltaX_fall;
                slope_fall.(mode{i})(r,2) = chan;
                slope_fall.(mode{i})(r,3) = stim;


                % Determine the rise time (N1 to P1)
                deltaY_rise = yP1-yN1;
                deltaX_rise = xP1 - xN1;
                slope_rise.(mode{i})(r,1) = deltaY_rise / deltaX_rise;
                slope_rise.(mode{i})(r,2) = chan;
                slope_rise.(mode{i})(r,3) = stim;

                

            end
          
          end
                    
        end
       
    end
    
end

% Check whether prop en clin have the same size
% Find zeros in (slope_fall.SPES_prop)          % Verwijderen
if sum(ismember(slope_fall.SPES_prop(:,1) ,0)) > 0
    rem_row = find(ismember(slope_fall.SPES_prop(:,1) ,0));
    
    slope_fall.SPES_clin(rem_row,:) = [];
    slope_fall.SPES_prop(rem_row,:) = [];
    slope_rise.SPES_clin(rem_row,:) = [];
    slope_rise.SPES_prop(rem_row,:) = [];
    
elseif sum(ismember(slope_rise.SPES_prop(:,1),0)) > 0
    rem_row = find(ismember(slope_rise.SPES_prop(:,1) ,0));
    
    slope_fall.SPES_clin(rem_row,:) = [];
    slope_fall.SPES_prop(rem_row,:) = [];
    slope_rise.SPES_clin(rem_row,:) = [];
    slope_rise.SPES_prop(rem_row,:) = [];

elseif sum(ismember(slope_fall.SPES_clin(:,1),0)) > 0
    rem_row = find(ismember(slope_fall.SPES_clin(:,1) ,0));
    
    slope_fall.SPES_clin(rem_row,:) = [];
    slope_fall.SPES_prop(rem_row,:) = [];
    slope_rise.SPES_clin(rem_row,:) = [];
    slope_rise.SPES_prop(rem_row,:) = [];
    
elseif sum(ismember(slope_rise.SPES_clin(:,1),0)) > 0
    rem_row = find(ismember(slope_rise.SPES_(:,1) ,0));
    
    slope_fall.SPES_clin(rem_row,:) = [];
    slope_fall.SPES_prop(rem_row,:) = [];
    slope_rise.SPES_clin(rem_row,:) = [];
    slope_rise.SPES_prop(rem_row,:) = [];
    
end

% Determine the median rise and fall times
N1_fall.median_fall_prop = median(slope_fall.SPES_prop(:,1));
N1_fall.fall_prop = slope_fall.SPES_prop(:,1);

N1_rise.median_rise_prop = median(slope_rise.SPES_prop(:,1));
N1_rise.rise_prop = slope_rise.SPES_prop(:,1);

N1_fall.median_fall_clin = median(slope_fall.SPES_clin(:,1));
N1_fall.fall_clin = slope_fall.SPES_clin(:,1);

N1_rise.median_rise_clin = median(slope_rise.SPES_clin(:,1));
N1_rise.rise_clin = slope_rise.SPES_clin(:,1);


% Save N1 rise and fall times for mathematical analysis
targetFolder = fullfile(myDataPath.CCEPpath, 'Visualise_agreement/N1_compare/N1_rise_and_fall/');
fileName_rise = [dataBase(1).sub_label,'_N1_rise.mat'];    
fileName_fall = [dataBase(1).sub_label,'_N1_fall.mat'];    

% Create the folder if it doesn't exist already.
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

save([targetFolder,fileName_rise], 'N1_rise');
save([targetFolder,fileName_fall], 'N1_fall');


end

%% Script to visualise the P0, N1 and P1

 % Find P0 and P1. Then the rise and fall times for the N1 peak
            % can be determined.
            
%             %%%%%%%% P0 %%%%%%%%%%%%%%%%%%%%
% % 
%             figure('Position',[391,61,1076,712])
%             plot(squeeze(dataBase_prop.cc_epoch_sorted_avg(chan,stim,:)),'k','LineWidth',2.5)
%             hold on
%             title(sprintf('SPES clin, %s, %s, %s',dataBase_clin.sub_label, dataBase_clin.stimpnames_avg{stim},  dataBase_clin.ch{chan}))
%             xlabel('time (samples)'); 
%             xlim([2*fs-10 (2*fs)+205]);
% %             
%             % Plot N1 detected by detector
%             plot(dataBase_prop.ccep.n1_peak_sample(chan,stim) ,dataBase_prop.ccep.n1_peak_amplitude(chan,stim)-2,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',15)
%             text(dataBase_clin.ccep.n1_peak_sample(chan,stim),dataBase_clin.ccep.n1_peak_amplitude(chan,stim)+100,'N1','FontSize',12)
%             hold on 
% %             ylabel('Potential \muV');
%             
%             % Create patch to indicate the 9 ms interval
%             ax = gca;
%             ylimits = ax.YTick;
%             patch([2*fs 2*fs+19 2*fs+19 2*fs],[min(ylimits) min(ylimits) max(ylimits) max(ylimits)],[0.6,0.2,0.2], 'EdgeAlpha',0); alpha(0.1) ;               
% %            
            % Find highest peak before N1 (P0)
            % The interval in wich de P0 is found is [9 ms (4115 samples) : % N1 location]
            % 4115 = 2*fs(2048) + 19 samples (stim artefact)
%             [yP0,xP0] = findpeaks(squeeze(dataBase_clin.cc_epoch_sorted_avg(chan,stim,[4115:dataBase_clin.ccep.n1_peak_sample(chan,stim)])));
%             
%             % When multiple peaks are found, select the highest (y)
%             if size(yP0,1) > 1    
%                [yP0,y0_loc] = (max(yP0) );
%                xP0 = xP0(y0_loc);
%             end
%             
%             % When no peak is found, take the point at 9 ms
%             if isempty(xP0)
%                 yP0 = dataBase_clin.cc_epoch_sorted_avg(chan,stim,2*fs+19);
%                 xP0 = 0;
%             end

%             % Plot the highest peak before N1 to check
%             hold on
%             plot((2*fs+19)+xP0,yP0,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
%             text((2*fs+19)+xP0,yP0+50,'P0','FontSize',12)

            
            %%%%%%%% P1 %%%%%%%%%%%%%%%%%%%%
             % Find highest peak before N1 (P1)
            % The interval in wich de P1 is found is [N1 location:100 ms (2*fs+205 samples = 4301)]
%             % 4115 = 2*fs(2048) + 19 samples (stim artefact)
%             [yP1,xP1] = findpeaks(squeeze(dataBase_clin.cc_epoch_sorted_avg(chan,stim,[dataBase_clin.ccep.n1_peak_sample(chan,stim):4301])));
%             
%             % When multiple peaks are found, select the first peak (x)
%             % or highest (y)???
%             if size(xP1,1) > 1    
%                [xP1,x1_loc] = (min(xP1) );
%                yP1 = yP1(x1_loc);
%             end

            % Plot the highest peak before N1 to check
%             hold on
%             plot((dataBase_clin.ccep.n1_peak_sample(chan,stim)+xP1),yP1,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15)
%             text((dataBase_clin.ccep.n1_peak_sample(chan,stim)+xP1),yP1-100,'P1','FontSize',12)
            

