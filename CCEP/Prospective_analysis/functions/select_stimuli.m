function remove_eventRow =  select_stimuli(dataBase, cfg,ll)

epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;
tt = (1:epoch_length*dataBase.ccep_header.Fs)/dataBase.ccep_header.Fs - epoch_prestim;

                %for ll = 1:size(dataBase.cc_stimsets_all,1)                % Find all the events in which this stimpair occurs

                    eventnum = find(strcmp(dataBase.tb_events.electrical_stimulation_site,[dataBase.cc_stimchans_all{ll,1}, '-',dataBase.cc_stimchans_all{ll,2}]));
                   
%                     sprintf('%s is stimulated %d time(s)', dataBase.stimpnames_all{ll},size(eventnum,1))
                    cc_epoch_sorted_all = [];
                    tt_epoch_sorted_all = [];
                    plot_ccep_prop =[];
                    stimnames_check= [];
                    figure('Position',[4,118,1914,938])
                    
                    for n = 1:size(eventnum,1)
                        % Propofol SPES
                        % Electrodes in the stimulation pair
                       
                        cc_epoch_sorted_all(:,:) = dataBase.data(:,dataBase.tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase.ccep_header.Fs)+1:dataBase.tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase.ccep_header.Fs));
                        tt_epoch_sorted_all(:,:) = dataBase.tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase.ccep_header.Fs)+1:dataBase.tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase.ccep_header.Fs);
               
                        % Remove electrodes which are stimulated
                        stim_elec_prop = dataBase.cc_stimsets_all(ll,:);
                        stimnames_check = dataBase.stimpnames_all(stim_elec_prop);
                       
                        plot_ccep_prop(:,:) = cc_epoch_sorted_all(:,:);
                        plot_ccep_prop(stim_elec_prop,:) = NaN;

                        % When the status of the channel is bad, not visualise the response
                        bad_sig(:,1) = find(contains(dataBase.tb_channels.status, {'bad'}));
                        plot_ccep_prop(bad_sig,:) = NaN;                       
                        plot_ccep_prop(:,tt>-0.01 & tt<0.02) = NaN;

                        % plot propofol SPES
                        subplot(1,size(eventnum,1),n)
                        hold on
                        plot(tt, plot_ccep_prop' + [0:1000:size(plot_ccep_prop,1)*1000-1])

                        str = sprintf('%s Propofol SPES', dataBase.stimpnames_all{ll});
                        title(str)
                        set(gca,'YTick',[0:1000:size(plot_ccep_prop,1)*1000-1],'YTickLabel',dataBase.ch) ;                  
                        xlim([-.2 0.3])
                        ylim([-1000 size(plot_ccep_prop,1)*1000-1])
                        ylabel('All stimuli of this stimulation pair' )
                        xlabel('time (s)')                           
                    end     
                        
                   remove_event = input('Which stimulus do you want to remove, type figurenumber like: [,], or press enter when you want to keep all: ')
                   remove_eventRow = eventnum(remove_event(:)) 


close all
                
end
