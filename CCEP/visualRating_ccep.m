
function dataBase = visualRating_ccep(dataBase,cfg)
fs = dataBase(1).ccep_header.Fs;

tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

dataBase.ccep.checked = zeros(size(dataBase.cc_epoch_sorted_avg,1),size(dataBase.cc_epoch_sorted_avg,2));

for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)
    for chan =1 :size(dataBase.cc_epoch_sorted_avg,1)
        
        if ~isnan(dataBase.ccep.n1_peak_sample(chan,stimp))
            % figure with left the epoch, and right zoomed in
            H=figure(1);
            H.Units = 'normalized';
            H.Position = [0.13 0.11 0.77 0.8];
            
            subplot(1,2,1)
            plot(tt,squeeze(dataBase.cc_epoch_sorted(chan,:,stimp,:)),'r');
            hold on
            plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            hold off
            xlim([-2 2])
            ylim([-2000 2000])
            xlabel('time(s)')
            ylabel('amplitude(uV)')
            title(sprintf('Electrode %s, stimulating %s-%s',dataBase.ch{chan},dataBase.cc_stimchans{stimp,1},dataBase.cc_stimchans{stimp,2}))
            
            subplot(1,2,2)
            plot(tt,squeeze(dataBase.cc_epoch_sorted(chan,:,stimp,:)),'r');
            hold on
            plot(tt,squeeze(dataBase.cc_epoch_sorted_avg(chan,stimp,:)),'k','linewidth',2);
            hold off
            xlim([-0.5 1])
            ylim([-750 750])
            title('Zoomed average signal')
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
            x = input('ER? [y/n] ','s');
            if strcmp(x,'y')
                dataBase.ccep.checked(chan,stimp) = 1 ;
            else
                dataBase.ccep.checked(chan,stimp) = 0 ;
            end
        end
        
    end
end