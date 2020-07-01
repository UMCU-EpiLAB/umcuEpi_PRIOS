function plot_ccep_av_stimp(dataBase,dataBase2, myDataPath, stimchans, LocOnes, TotOnesStim, dif_mat)

LocaOnes = LocOnes{:,:};             % Create matrix of table

[ER_in10(:,2), ER_in10(:,1)] = find(dif_mat == 1);              % [stimpairs electrodes]

for i = 1:size(ER_in10(:,2))                                    % For the number of ones detected
    ER_in10st(i,1) = stimchans(ER_in10(i,1))';
    ER_in10st(i,2) = dataBase.ch(ER_in10(i,2));    
end


tt = dataBase.tt;    
   for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)            % for all stimulation pairs
        
         Stimpnm = stimchans{stimp};    
         counter = 0;
         counter_10 = 0;
         nmbr_of_ones = TotOnesStim(stimp);
         name = cell(1,nmbr_of_ones);
         ccep_plot = zeros(length(tt),1);
         ccep_plot2 = zeros(length(tt),1);
         NameER = cell(1,nmbr_of_ones);
              
         
        for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)            % for all electrodes
             elecnm = dataBase.ch{elec};
             
            for i = 1:size(LocaOnes)     
                if ismember({Stimpnm}, LocaOnes{i,1}) && ismember(dataBase.ch{elec}, LocaOnes(i,2)) 
                    counter = counter+1;
                    name(:,counter) =  [LocaOnes(i,2)];                 % names of electrodes with a 1 on stimulation pair stimp.
                     
                    ccep_plot(:,counter) = squeeze(dataBase.cc_epoch_sorted_avg(elec,stimp,:));
                    ccep_plot2(:,counter) = squeeze(dataBase2.cc_epoch_sorted_avg(elec,stimp,:));       % stimulations of the 2 stims
                    ccep_plot(tt>-0.010 & tt<0.010) = NaN;
                    ccep_plot2(tt>-0.010 & tt<0.010) = NaN;
                          
                    figure('Position',[1400 0 700 700])
                    %hold on
                    xlim([-.2 1.5])
                    ylabel('Average per electrodes (mV)')
                    xlabel('time (s)') 
                       
                    subplot(2,1,1);
                    plot(tt,  ccep_plot + [500:500:counter*500]);   
                    str = sprintf('Stimulation pair %s for 10 stimps', stimchans{stimp,1});
                    title(str)
                    xlim([-.2 1.5])
                    ylabel('Average per electrodes (mV)')
                    xlabel('time (s)') 
                    set(gca,'YTick',500*(1:counter),'YTickLabel',name(:)) ;
                    %hold  on

                    subplot(2,1,2);
                    plot(tt,  ccep_plot2 + [500:500:counter*500]);   
                    str2 = sprintf('Stimulation pair %s for 2 stimps', stimchans{stimp,1});
                    title(str2)
                    xlim([-.2 1.5])
                    ylabel('Average per electrodes (mV)')
                    xlabel('time (s)') 
                    set(gca,'YTick',500*(1:counter),'YTickLabel',name(:)) ;
                    %hold off
                
                
                 for j = 1:length(ER_in10st)
                    if ismember({Stimpnm}, ER_in10st{j,1}) && ismember(dataBase.ch{elec}, ER_in10st(j,2))   
                        counter_10 = counter_10+1;
                        NameER(:,counter_10) =  [ER_in10st(j,2)];
                        str_main = sprintf('10 stims evoke ER in %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s',NameER{:});        
                        sgtitle(str_main)
                    end
                 end
                

    % Save the figures
        if dataBase.save_fig==1
            % create folder to save figures
            if ~ exist(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label),'dir')

                mkdir(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label));
            end

            % filename
            figureName = fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,...
                [dataBase.sub_label '_stimp' stimchans{stimp}]);
            set(gcf,'PaperPositionMode','auto');
            print('-dpng','-r300',figureName);
        else
            pause
        end
          end
            
            
    close all
         

            end
     end              
    end
end

          
