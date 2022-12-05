function visualRating_ccep(dataBase, cfg, endstimp, myDataPath)
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

if any(contains(fieldnames(dataBase.ccep),'obs_tab') )         % if the N1-check has been performed before
    n1_peak_amplitude_check = dataBase.ccep.n1_peak_amplitude_check;
    n1_peak_sample_check = dataBase.ccep.n1_peak_sample_check;
    obs_tab = dataBase.ccep.obs_tab;         % Table is filled with the marker of the observer.

else
    n1_peak_amplitude_check = NaN(size(n1_peak_amplitude));
    n1_peak_sample_check = NaN(size(n1_peak_sample));
    obs_tab = cell(size(n1_peak_amplitude));         % Table is filled with the marker of the observer.
end

%% Check the automatically detected ER for every stimulation-electrode
% combination in which an N1 is detected.

n = numel(n1_peak_sample_check(:,1:endstimp))+1;

% When visual N1-check was started already, then start at the last stim
% checked (endstimp) +1
for stimp = endstimp+1:size(dataBase.cc_epoch_sorted_reref_avg,2)

    for chan = 1:size(dataBase.cc_epoch_sorted_reref_avg,1)

        if ~isnan(n1_peak_sample(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H = figure(1);
            H.Units = 'normalized';
            H.Position = [0.14,0.0625,0.77,0.7];

            subplot(1,2,1)
            plot(tt,squeeze(dataBase.cc_epoch_sorted_reref(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
            hold on
            plot(tt,squeeze(dataBase.cc_epoch_sorted_reref_avg(chan,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
            plot(tt(n1_peak_sample(chan,stimp)), n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
            hold off
            xlim([-1 1])
            ylim([-2000 1500])
            xlabel('Time(s)')
            ylabel('Potential (\muV)')
            title(sprintf('Electrode %s, stimulating %s-%s',dataBase.ch{chan},dataBase.cc_stimchans{stimp,:}))

            % Zoom in of the left figure
            subplot(1,2,2)
            p1 = plot(tt,squeeze(dataBase.cc_epoch_sorted_reref(chan,stimp,:,:)),'r:');                % plot the 10 separate stimulations
            hold on
            p2 = plot(tt,squeeze(dataBase.cc_epoch_sorted_reref_avg(chan,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
            p3 = plot(tt(n1_peak_sample(chan,stimp)), n1_peak_amplitude(chan,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4);
            hold off
            xlim([-0.05 0.2])

            % make sure the minimal value of the n1_peak_amplitude is
            % visible in the zoomed signal
            if n1_peak_amplitude(chan,stimp) > -1000
                ylim([-1000 450])
            else
                ymin = floor(1.1*n1_peak_amplitude(chan,stimp));
                ylim([ymin 450])
            end

            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Potential \muV')

            % Create a patch for the -1/5:9 ms interval in which no
            % physiological activity can be measured.
            patch([0 0.009 0.009 0],[-1200 -1200 1200 1200],[0.6,0.2,0.2],'EdgeAlpha',0)
            alpha(0.2)

            legend([p1(1),p2(1),p3(1)],'indiv responses','average','detected CCEP')

            % percentage of pictures that you have scored so far
            perc = n/size(n1_peak_amplitude(:),1)*100;

            fprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this an N1 (y/n)? [y/n] \n  If incorrect N1 detection, select correct N1 peak in right figure (zoomed) and press ENTER. \n',...
                perc,dataBase.cc_stimchans{stimp,:},dataBase.ch{chan});
            currkey = 0;
            cp =[];

            % select new N1 or categorize as good N1 or no N1
            % When incorrect N1 is selected, click on correct N1, a blue
            % spot will appear, then press enter. The new coordinates will
            % be entered in n1_peak_amplitude and sample.

            while ~strcmp(currkey,{'y','n','d',char(13)}) % y,n,d or enter

                w = waitforbuttonpress; % 0 = mouse, other = key
                if w == 0 % = mouse
                    % draw correct N1
                    cp = get(gca,'CurrentPoint');

                    % find sample number closest to the selected point
                    [~,sampnum] = min(abs(tt-cp(1,1)));

                    % find nearby peak (within 20ms of selected point)
                    [~,locs] = findpeaks(-1*squeeze(dataBase.cc_epoch_sorted_reref_avg(chan,stimp,...
                        sampnum-round(0.01*fs):sampnum+round(0.01*fs))),'NPeaks',1,'SortStr','descend');

                    locsamp = sampnum-round(0.01*fs)+locs-1;  % find x-position of nearby peak

                    hold on
                    plot(tt(locsamp),dataBase.cc_epoch_sorted_reref_avg(chan,stimp,locsamp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6); drawnow;
                    n1_peak_amplitude_check(chan,stimp) = dataBase.cc_epoch_sorted_reref_avg(chan,stimp,locsamp);
                    n1_peak_sample_check(chan,stimp) = locsamp ;
                    fprintf('sample: %1.4f, amplitude: %2.1f \n',tt(locsamp),dataBase.cc_epoch_sorted_reref_avg(chan,stimp,locsamp))

                    hold off

                    %%% Add marking of observer to a table
                    obs_tab{chan,stimp} = 'c';

                elseif w == 1 % keyboard
                    currkey = get(gcf,'CurrentCharacter');

                    if strcmp(currkey,'y') && isempty(cp) % if nothing is annotated,
                        n1_peak_amplitude_check(chan,stimp) = n1_peak_amplitude(chan,stimp) ;
                        n1_peak_sample_check(chan,stimp) = n1_peak_sample(chan,stimp) ;
                        
                        hold off
                    
                    elseif strcmp(currkey,'y') && ~isempty(cp) % if something is annotated
                        % do nothing because it is already saved correctly
                        % in n1_peak_sample_check and
                        % n1_peak_amplitude_check 

                    elseif (strcmp(currkey,'y') || strcmp(currkey,'d')) && ~isempty(cp) % if something is annotated
                        % do nothing because it is already saved correctly
                        % in n1_peak_sample_check and
                        % n1_peak_amplitude_check

                    elseif strcmp(currkey,'n')
                        n1_peak_amplitude_check(chan,stimp) = NaN ;
                        n1_peak_sample_check(chan,stimp) = NaN ;
                        
                        hold off

                        currkey = 'n'; % enter is interpreted as 'n', and displayed as 'n' in line obs_tab

                    end

                    % Add marking of observer to a table
                    obs_tab{chan,stimp} = currkey;

                    if ~strcmp(currkey,{'y','n',char(13)})
                        fprintf('---- ANSWER ---- : %s, Select one of the options [y/n/enter]: \n',currkey);

                    elseif strcmp(currkey,{char(13)})
                        fprintf('---- ANSWER ---- : new N1 was selected: \n');

                    elseif strcmp(currkey,{'y','n'})
                        fprintf('---- ANSWER ---- : %s \n \n',currkey);

                    end % if there was no y, n, space entered on the keyboard
                end % if there was a mouse click
            end % while no y,n,d, or space was pressed on the keyboard
        end % if N1 was detected
        
        hold off

        n = n+1;

    end % for-loop channel

    % save also till which stimpair visual N1s are checked.
    ccep.checkUntilStimp = stimp;
    ccep.n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
    ccep.n1_peak_sample = dataBase.ccep.n1_peak_sample;
    ccep.n1_peak_amplitude_check = n1_peak_amplitude_check;
    ccep.n1_peak_sample_check = n1_peak_sample_check;
    ccep.obs_tab = obs_tab;

    ccep.amplitude_thresh = dataBase.ccep.amplitude_thresh;
    ccep.n1_peak_range = dataBase.ccep.n1_peak_range;
    ccep.cc_stimchans = dataBase.cc_stimchans;
    ccep.cc_stimsets = dataBase.cc_stimsets;
    ccep.dataName = dataBase.dataName;
    ccep.ch = dataBase.ch;
    ccep.tt = dataBase.tt;
    ccep.epoch_length = cfg.epoch_length;
    ccep.epoch_prestim = cfg.epoch_prestim;

    filename = [dataBase.sub_label,'_',dataBase.ses_label,'_', ...
        dataBase.task_label,'_N1sChecked_', cfg.observer, '.mat'];

    filefolder = fullfile(myDataPath.CCEPpath,'checkedN1s');

    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end

    % save file during scoring in case of error
    save(fullfile(filefolder,filename),'ccep','-struct');
    
end % for-loop stimulus pair

fprintf('CCEPs are saved for %s for subject %s \n' , dataBase(1).task_label, dataBase(1).sub_label);

end

