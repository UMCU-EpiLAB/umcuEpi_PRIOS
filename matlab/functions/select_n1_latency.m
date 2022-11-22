% Observers are able to select another N1-peak than found by the detector. 
% When two observers have selected two different points, take the mean of the two observed points when the difference is <
% 0.03 seconds. When the difference is more, than perform an visual check and
% select the right point. 
            
function dataBase = select_n1_latency(rater1, rater2, myDataPath, sublabel, tasklabel, subj)
      
% load CCEP's to plot later
files_in_folder = dir([myDataPath.dataPath,'derivatives','/CCEPs']);
file_name = [['sub-',sublabel],'_',['task-',tasklabel],'_CCEP.mat'];

for i = 1:size(files_in_folder,1)
    if contains(files_in_folder(i).name,file_name)
        
        if isequal(tasklabel , 'SPESclin')
            tic
            load([[myDataPath.dataPath,'derivatives','/CCEPs/'],files_in_folder(i).name]); % this could take more than 2 minutes, only do ones and save.
            toc
            ccep = dataBase_clin;
        elseif  isequal(tasklabel , 'SPESprop')
            tic
            load([[myDataPath.dataPath,'derivatives','/CCEPs/'],files_in_folder(i).name]); % this could take more than 2 minutes, only do ones and save.
            toc
            ccep = dataBase_prop;
        end

       
    end

end
% preallocate
mean_N1_lat = NaN(size(ccep.ch,1), size(ccep.cc_stimsets_avg,1));

for stimp = 1:size(rater1.ccep.n1_peak_sample_check,2)

    for elec = 1:size(rater1.ccep.n1_peak_sample_check,1)
        if ~isnan(rater1.ccep.n1_peak_sample_check(elec,stimp)) && ~isnan(rater2.ccep.n1_peak_sample_check(elec,stimp))
            % rater 1 and rater 2 should have the same nan and non-nan
            % because in script interobserver_analysis only responses with
            % an N1 scored by both observers are saved, when only 1
            % observer scored an N1 then it is replaced with NaN for both
            % observers.

            lat_R1 = rater1.ccep.n1_peak_sample_check(elec,stimp);
            lat_R2 = rater2.ccep.n1_peak_sample_check(elec,stimp);

            if abs(lat_R1 - lat_R2) > 5 
                % plot figure to determine whether observer 1 or 2 was right
                H=figure(1);
                H.Units = 'normalized';
                H.Position = [0.14,0.0625,0.77,0.7];
    
                plot(ccep.tt, squeeze(ccep.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
                hold on
                
                plot(ccep.tt(rater1.ccep.n1_peak_sample_check(elec,stimp)), rater1.ccep.n1_peak_amplitude_check(elec,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
                plot(ccep.tt(rater2.ccep.n1_peak_sample_check(elec,stimp)), rater2.ccep.n1_peak_amplitude_check(elec,stimp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',4)
                plot(ccep.tt,squeeze(ccep.cc_epoch_sorted_select_reref(elec,stimp,:,:)),'r:');                % plot the 10 separate stimulations

                % Create a patch for the -1/5:9 ms interval in which no
                % physiological activity can be measured.
                patch([0 0.009 0.009 0],[-800 -800 750 750],[0.6,0.2,0.2],'EdgeAlpha',0)
                alpha(0.2)
                
                hold off
    
                xlim([-0.05 0.2])
                ymin = floor(1.1*rater1.ccep.n1_peak_amplitude(elec,stimp));
                ylim([ymin 600])
                xlabel('Time(s)')
                ylabel('Potential (\muV)')
                title(sprintf('Electrode %s, stimulating %s',ccep.ch{elec},ccep.stimpnames_avg{stimp}))
                legend(' ','Rater1','Rater2',' ')


                % Decide which of the two N1-peaks is the right one
                prompt = 'Which N1-peak is correct selected? Rater1(red) or Rater2(blue) 1/2: ';
                str = input(prompt,'s');

                if isequal(str,'1')
                    mean_N1_lat(elec,stimp) = rater1.ccep.n1_peak_sample_check(elec,stimp);

                elseif isequal(str,'2')
                    mean_N1_lat(elec,stimp) = rater2.ccep.n1_peak_sample_check(elec,stimp);

                end

            else % when latencies for both observers was (almost) equal  
                mean_N1_lat(elec,stimp) = round(mean([rater1.ccep.n1_peak_sample_check(elec,stimp),rater2.ccep.n1_peak_sample_check(elec,stimp)]));

            end

        end

    end

end


    
if isequal(tasklabel,'SPESclin')  
  
    dataBase.ccep_clin.n1_peak_sample = mean_N1_lat;        
    dataBase.ccep_clin.sub_label = sublabel;
    dataBase.ccep_clin.task_label= tasklabel;
    dataBase.ccep_clin.stimsets_avg = rater1.ccep.stimsets_avg;
    dataBase.ccep_clin.stimpnames_avg = rater1.ccep.stimpnames_avg;
    dataBase.ccep_clin.ch = rater1.ccep.ch;

    % Save file in inter_observer folder REPLACE
    % Keep the original in the original patient sub folders.
    filefolder = fullfile(myDataPath.CCEP_interObVar);
    filename = [['sub-',dataBase.ccep_clin.sub_label],'_ses-1_',['task-',dataBase.ccep_clin.task_label],'_meanN1Lat.mat'];
    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end

    % Save to make sure that you do not have to select the right N1-peaks
    % when two observers selected different peaks
    save(fullfile(filefolder,filename),'dataBase');

elseif isequal(tasklabel,'SPESprop')
    
    dataBase.ccep_prop.n1_peak_sample = mean_N1_lat;        
    dataBase.ccep_prop.sub_label = sublabel;
    dataBase.ccep_prop.task_label= tasklabel;
    dataBase.ccep_prop.stimsets_avg = rater1.ccep.stimsets_avg;
    dataBase.ccep_prop.stimpnames_avg = rater1.ccep.stimpnames_avg;
    dataBase.ccep_prop.ch = rater1.ccep.ch;


    % Save file in inter_observer folder REPLACE
    % Keep the original in the original patient sub folders.
    filefolder = fullfile(myDataPath.CCEP_interObVar);
    filename = [['sub-',dataBase.ccep_prop.sub_label],'_ses-1_',['task-',dataBase.ccep_prop.task_label],'_meanN1Lat.mat'];
    if ~exist(filefolder,'dir')
        mkdir(filefolder)
    end
    
    % Save to make sure that you do not have to select the right N1-peaks
    % when two observers selected different peaks
    save(fullfile(filefolder,filename),'dataBase');


end



end



