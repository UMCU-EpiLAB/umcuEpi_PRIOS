function plot_ccep_av_stimp2(dataBase,dataBase2, myDataPath, stimchans, LocOnes, TotOnesStim, dif_mat)

LocaOnes = LocOnes{:,:};             % Create matrix of table

[ER_in10(:,2), ER_in10(:,1)] = find(dif_mat == 1);              % [stimpairs electrodes]

for i = 1:size(ER_in10(:,2))                                    % For the number of ones detected
    ER_in10st(i,1) = stimchans(ER_in10(i,1))';
    ER_in10st(i,2) = dataBase.ch(ER_in10(i,2));    
end


tt = dataBase.tt;    
   for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)            % for all stimulation pairs
        
        Stimpnm = stimchans{stimp};
        figure('Position',[1400 0 700 700])
        hold on
        xlim([-.2 1.5])
        ylabel('Average per electrodes (mV)')
        xlabel('time (s)') 
                
         counter = 0;
         nmbr_of_ones = TotOnesStim(stimp);
         name = cell(1,nmbr_of_ones);
            
       
        for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)            % for all electrodes
             elecnm = dataBase.ch{elec};
             
            for i = 1:size(LocaOnes)     
                if ismember({Stimpnm}, LocaOnes{i,1}) && ismember(dataBase.ch{elec}, LocaOnes(i,2)) 
                     counter = counter+1;
                     name(:,counter) =  [LocaOnes(i,2)];                 % names of electrodes with a 1 on stimulation pair stimp.
%                 end
%             end
%         end
%         
%         
         for nmr_plot = 1:nmbr_of_ones                            % total number of ones for stimulation pair stimp  
%             Electodes = [dataBase.ch(:)];
%         
%             loc_electr = all(ismember(Electodes, name{nmr_plot}),2 ) ; 
%             loc_elect(:,nmr_plot) = find(loc_electr) ;              % Number of the location of the electrode with a 1 on stimpair stimp.
         
                ccep_plot = squeeze(dataBase.cc_epoch_sorted_avg(elec,stimp,:));
                ccep_plot2 = squeeze(dataBase2.cc_epoch_sorted_avg(elec,stimp,:));       % stimulations of the 2 stims
                ccep_plot(tt>-0.010 & tt<0.010) = NaN;
                ccep_plot2(tt>-0.010 & tt<0.010) = NaN;

                subplot(2,1,1);
                plot(tt, nmr_plot*500+ccep_plot, 'LineWidth', 2);             % plot(tt, loc_elect*500+ccep_plot, 'k', 'LineWidth', 3);  
                str = sprintf('Stimulation pair %s for 10 stimps', stimchans{stimp,1});
                title(str)
                xlim([-.2 1.5])
                ylabel('Average per electrodes (mV)')
                xlabel('time (s)') 
                set(gca,'YTick',500*(1:counter),'YTickLabel',name(:)) ;
                hold on 
                
                
                subplot(2,1,2);
                plot(tt, nmr_plot*500+ccep_plot2, 'LineWidth', 2); 
                str2 = sprintf('Stimulation pair %s for 2 stimps', stimchans{stimp,1});
                title(str2)
                xlim([-.2 1.5])
                ylabel('Average per electrodes (mV)')
                xlabel('time (s)') 
                set(gca,'YTick',500*(1:counter),'YTickLabel',name(:)) ;
                hold on

                 for j = 1:length(ER_in10st)
                    if ismember({Stimpnm}, ER_in10st{j,1}) && ismember(dataBase.ch{elec}, ER_in10st(j,2))   
                        %% !!!! Naam wordt nu telkens overschreven als er meerdere zijn ER in 10 stims
                        str_main = sprintf('ER in 10 stims. Stimulation pair %s for electrode %s', Stimpnm, elecnm);        
                        sgtitle(str_main)
                    end
                 end
         end
           
        end
       
     hold off

    % Save the figures
        if dataBase.save_fig==1
            % create folder to save figures
            if ~ exist(fullfile(myDataPath.CCEPpath,'av_ccep_figures2',dataBase.sub_label),'dir')

                mkdir(fullfile(myDataPath.CCEPpath,'av_ccep_figures2',dataBase.sub_label));
            end

            % filename
            figureName = fullfile(myDataPath.CCEPpath,'av_ccep_figures2',dataBase.sub_label,...
                [dataBase.sub_label '_' dataBase.NmbrofStims '_stimp' stimchans{stimp}]);
            set(gcf,'PaperPositionMode','auto');
            print('-dpng','-r300',figureName);
        else
            pause
        end
    close all
   
            end        
    end
end

          
