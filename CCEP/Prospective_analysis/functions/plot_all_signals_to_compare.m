
 %%% Propofol SPES
 stimp = 40;        % PRIOS06 C15C16 = 39, 40
 figure()
 
bad_sig_prop(:,1) = find(contains(dataBase_prop.tb_channels.status, {'bad'}));
dataBase_prop.cc_epoch_sorted(bad_sig_prop,:,:) = NaN;

titleRows = ['geen suppre';'suppression';'tussen supp';'weer normal'];  

 for turn = 1:4
        % Electrodes in the stimulation pair
        stim_elec_prop = dataBase_prop.cc_stimsets_all(stimp,:);
        
        % Remove electrodes which are stimulated
        plot_ccep_prop = dataBase_prop.cc_epoch_sorted(:,:,:,:);
        plot_ccep_prop(stim_elec_prop,:,:) = NaN;

        plot_ccep_prop = squeeze(plot_ccep_prop(:,turn,stimp,:));
        plot_ccep_prop(:,tt>-0.01 & tt<0.02) = NaN;
        
     
        % plot propofol SPES
        subplot(1,4,turn)
        plot(tt, plot_ccep_prop' + [0:1000:size(plot_ccep_prop,1)*1000-1])
        
        set(gca,'YTick',[0:1000:size(plot_ccep_prop,1)*1000-1],'YTickLabel',dataBase_prop.ch) ;                  
        xlim([-.2 1.5])
        ylim([-1000 size(plot_ccep_prop,1)*1000-1])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 
        title(sprintf('C16C15 %s',titleRows(turn,:)))

 end
 
