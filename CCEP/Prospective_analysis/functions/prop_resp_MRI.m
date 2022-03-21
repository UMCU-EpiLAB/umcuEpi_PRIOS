function  prop_resp_MRI(dataBase, myDataPath)

% Observers are able to select another N1-peak than found by the detector. 
% When two observers have selected two different points, take the mean of the two observed points when the difference is <
% 0.03 seconds. When the difference is more, than perform an visual check and
% select the right point. 
          

 for pat = 1:size(dataBase,2)   
       
         % load CCEP's to plot later
        files_in_folder = dir([myDataPath.dataPath,'derivatives','/CCEPs']);
        file_name_clin = [['sub-',dataBase(pat).ccep_clin.sub_label  ],'_',['task-',dataBase(pat).ccep_clin.task_label],'_CCEP.mat'];
        file_name_prop = [['sub-',dataBase(pat).ccep_clin.sub_label  ],'_',['task-',dataBase(pat).ccep_prop.task_label],'_CCEP.mat'];


        for i = 1:size(files_in_folder,1)
            if contains(files_in_folder(i).name,file_name_clin)

                    load([[myDataPath.dataPath,'derivatives','/CCEPs/'],files_in_folder(i).name]); % this could take more than 2 minutes, only do ones and save.
                
            elseif  contains(files_in_folder(i).name,file_name_prop)
                   
                    load([[myDataPath.dataPath,'derivatives','/CCEPs/'],files_in_folder(i).name]); % this could take more than 2 minutes, only do ones and save.

            end
        
        end
        
        
        for p = 1:size(dataBase(pat).elec_in_prop,1)

            stimp = dataBase(pat).elec_in_prop(p,1);
            elec = dataBase(pat).elec_in_prop(p,2);

            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.14,0.0625,0.77,0.7];

            %plot clinical SPES left
            subplot(1,2,1)
            plot(dataBase_clin.tt, squeeze(dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
            hold on
            plot(dataBase_clin.tt,squeeze(dataBase_clin.cc_epoch_sorted_select_reref(elec,stimp,:,:)),'r:');                % plot the 10 separate stimulations
           
            if ~isnan(dataBase(pat).ccep_clin.n1_sample_ori(elec,stimp))
                plot(dataBase_prop.tt(dataBase(pat).ccep_clin.n1_sample_ori(elec,stimp)), dataBase(pat).ccep_clin.n1_amp_ori(elec,stimp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',4)
            end

            hold off

            xlim([-0.05 0.2])
            ylim([-600 600])
            xlabel('Time(s)')
            ylabel('Potential (\muV)')
            title(sprintf('CLIN - Electrode %s, stimulating %s',dataBase_clin.ch{elec},dataBase_clin.stimpnames_avg{stimp}))

            subplot(1,2,2)
            plot(dataBase_prop.tt, squeeze(dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
            hold on
            plot(dataBase_prop.tt,squeeze(dataBase_prop.cc_epoch_sorted_select_reref(elec,stimp,:,:)),'r:');                % plot the 10 separate stimulations
            plot(dataBase_prop.tt(dataBase_prop.ccep.n1_peak_sample(elec,stimp)), dataBase_prop.ccep.n1_peak_amplitude(elec,stimp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',4)


            xlim([-0.05 0.2])
            ylim([-600 600])
            xlabel('Time(s)')
            ylabel('Potential (\muV)')
            title(sprintf('PROP - Electrode %s, stimulating %s',dataBase_clin.ch{elec},dataBase_clin.stimpnames_avg{stimp}))              
            

            filefolder = fullfile(myDataPath.CCEPpath,'Response_in_prop/');
            if ~isnan(dataBase(pat).ccep_clin.n1_sample_ori(elec,stimp))
                fileName = sprintf('%s,%s,%s',dataBase_clin.sub_label ,dataBase_clin.ch{elec},dataBase_clin.stimpnames_avg{stimp});
            else
                fileName = sprintf('CLIN not detected by detector %s,%s,%s',dataBase_clin.sub_label ,dataBase_clin.ch{elec},dataBase_clin.stimpnames_avg{stimp});

            end

            
            % Save to make sure that you do not have to select the right N1-peaks
            % when two observers selected different peaks
            set(gcf,'PaperPositionMode','auto')
            print('-dpng','-r300',[filefolder,fileName])

            hold off;
        end
 end

