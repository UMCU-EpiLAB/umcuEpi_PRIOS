function dataBase = visualRating_ccep(dataBase, cfg, endstimp, myDataPath)
% INSTRUCTIONS
% select new point: select new point > enter or y
% correct: y
% incorrect: n or enter

close all

tt = dataBase.tt;
fs = dataBase.ccep_header.Fs;  

% Load the N1 amplitude and sample number of automatic detection
n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
n1_peak_sample = dataBase.ccep.n1_peak_sample;


if sum(strcmp(fieldnames(dataBase.ccep),'check')) == 1          % if the N1-check has been performed before
    ccep.n1_peak_amplitude_check = dataBase.ccep.n1_peak_amplitude_check;
    ccep.n1_peak_sample_check = dataBase.ccep.n1_peak_sample_check;
else
    ccep.n1_peak_amplitude_check = NaN(size(n1_peak_amplitude));
    ccep.n1_peak_sample_check = NaN(size(n1_peak_sample));
end



%% Check the automatically detected ER for every stimulation-electrode
% combination in which an N1 is detected.

% Preallocation
obs_tab = cell(size(n1_peak_amplitude,1),size(n1_peak_amplitude,2));         % Table is filled with the marker of the observer. to save doubtfull obvervations


% When visual N1-check was started already, then start at the last stim
% checked (endstimp) +1
for stimp = endstimp+1:size(dataBase.cc_epoch_sorted_avg,2)
        
    for chan =1 :size(dataBase.cc_epoch_sorted_avg,1)
        
        if ~isnan(dataBase.ccep.n1_peak_sample(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.31 0.77 0.7];
            
            subplot(1,2,1)            
            if strcmp(cfg.reref,'y')
                plot(tt,squeeze(dataBase.cc_epoch_sorted_select_reref(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
                hold on
                plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k-.','linewidth',1);      % plot the not-rereference signal in a dashed line
                plot(tt,squeeze(dataBase.cc_epoch_sorted_select_reref_avg(chan,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
                
            else
                plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
                hold on
                plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
                
            end
            
            plot(tt(n1_peak_sample(chan,stimp)), n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
            hold off
            xlim([-1 1])
            xlabel('Time(s)')
            ylabel('Potential (\muV)')
            title(sprintf('Electrode %s, stimulating %s',dataBase.ch{chan},dataBase.stimpnames_avg{stimp}))
            
                        
            % Zoom in of the left figure
            subplot(1,2,2)
            if strcmp(cfg.reref,'y')
                p1 = plot(tt,squeeze(dataBase.cc_epoch_sorted_select_reref(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
                hold on
                p2 = plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k-.','linewidth',1);      % plot the not-rereference signal in a dashed line
                p3 = plot(tt,squeeze(dataBase.cc_epoch_sorted_select_reref_avg(chan,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line

            else
                p1 = plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
                p2 = plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
                
            end
         
            p4 = plot(tt(n1_peak_sample(chan,stimp)), n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4);
            hold off
            xlim([-0.05 0.2])
            ylim([-800 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Potential \muV')
            
            % Create a patch for the -1/5:9 ms interval in which no
            % physiological activity can be measured.
            patch([0 0.009 0.009 0],[-800 -800 750 750],[0.6,0.2,0.2],'EdgeAlpha',0)
            alpha(0.2)
            
            if strcmp(cfg.reref,'y')
                legend([p1(1),p2(1),p3(1),p4(1)],'indiv responses','average','average reref','detected CCEP')    
            elseif strcmp(cfg.reref,'n')
                legend([p1(1),p2(1),p4(1)],'indiv responses','average','detected CCEP')
            end
            
            currkey = 0;
            fprintf('N1 [y/n], if incorrect N1, select correct N1 and press enter. Only enter is equal to not correct \n')
            
            % select new N1 or categorize as good N1 or no N1
            % When incorrect N1 is selected, click on correct N1, a blue
            % stip will occure, then press enter! The new coordinates will
            % show in n1_peak_amplitude and sample.
            while ~strcmp(currkey,{'y','n',char(13)})
                cp =[];
                w = waitforbuttonpress; % 0 = mouse, other = key
                if w == 0
                    % draw correct N1
                    cp = get(gca,'CurrentPoint');
                                                              
                    % find sample number closest to the selected point
                    [~,sampnum] = min(abs(tt-cp(1,1)));
                    
                    % find nearby peak
                    if strcmp(cfg.reref,'y')
                        [~,locs] = findpeaks(-1*squeeze(dataBase.cc_epoch_sorted_select_reref_avg(chan,stimp,...
                            sampnum-round(0.01*fs):sampnum+round(0.01*fs))),'NPeaks',1,'SortStr','descend');
                    else
                        [~,locs] = findpeaks(-1*squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp, sampnum-round(0.01*fs):sampnum+round(0.01*fs))),...
                            'NPeaks',1,'SortStr','descend');                        
                    end
                    
                    
                    % find x-position of nearby peak
                    locsamp = sampnum-round(0.01*fs)+locs-1;
                    
                    hold on
                    if strcmp(cfg.reref,'y')
                        plot(tt(locsamp),dataBase.cc_epoch_sorted_select_reref_avg(chan,stimp,locsamp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6); drawnow;

                    else
                        plot(tt(locsamp),dataBase.cc_epoch_sorted_avg(chan,stimp,locsamp),'bo','MarkerFaceColor','b','MarkerSize',4); drawnow;

                    end
                         
                    ccep.n1_peak_sample_check(chan,stimp) = locsamp ;
                    
                    if strcmp(cfg.reref,'y')
                         ccep.n1_peak_amplitude_check(chan,stimp) = dataBase.cc_epoch_sorted_select_reref_avg(chan,stimp,locsamp);
                    else
                         ccep.n1_peak_amplitude_check(chan,stimp) = dataBase.cc_epoch_sorted_avg(chan,stimp,locsamp) ;
                    end
                    
                   hold off
                   
                   %%% Add marking of observer to a table
                   obs_tab{chan,stimp} = 'c';
                   
                   
                elseif w == 1
                    currkey = get(gcf,'CurrentCharacter');
                    
                    if strcmp(currkey,'y') || strcmp(currkey,'d') && isempty(cp)
                        ccep.n1_peak_amplitude_check(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        ccep.n1_peak_sample_check(chan,stimp) = n1_peak_sample(chan,stimp) ;
                        hold off
                   
                    elseif strcmp(currkey,'n')              
                        ccep.n1_peak_amplitude_check(chan,stimp) = NaN ;
                        ccep.n1_peak_sample_check(chan,stimp) = NaN ;
                        hold off
                    end
                    
                     %%% Add marking of observer to a table
                    obs_tab{chan,stimp} = currkey; 
                        
                end
            end
            
        end
       hold off
    end
    
   % save also till which stimpair visual N1s are checked.
    ccep.checkUntilStimp = stimp;
    ccep.n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
    ccep.n1_peak_sample = dataBase.ccep.n1_peak_sample;
    
    filename = [dataBase.sub_label,'_',dataBase.ses_label,'_',dataBase.task_name,'_N1sChecked.mat'];
    filefolder = fullfile(myDataPath.CCEPpath, dataBase.sub_label, dataBase.ses_label, dataBase.task_name,'/');
    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end

    % save file during scoring in case of error
    save(fullfile(filefolder,filename),'-struct','ccep');
end

dataBase.ccep.n1_peak_amplitude_check = ccep.n1_peak_amplitude_check;
dataBase.ccep.n1_peak_sample_check = ccep.n1_peak_sample_check;
dataBase.ccep.n1_peak_amplitude = n1_peak_amplitude;
dataBase.ccep.n1_peak_sample = n1_peak_sample;
dataBase.ccep.checkUntilStimp = stimp;


end

