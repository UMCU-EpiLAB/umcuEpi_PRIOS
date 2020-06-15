function plot_ccep_av_stimp(dataBase, myDataPath, stimchans, LocOnes, TotOnesStim)

LocaOnes = LocOnes{:,:};             % Create matrix of table

tt = dataBase.tt;    
   for stimp = 1:size(dataBase.cc_epoch_sorted_avg,2)            % for all stimulation pairs
        
        Stimpnm = stimchans{stimp};
        figure('Position',[200 0 700 700])
        hold on

        xlim([-.2 1.5])
        ylabel('Average per electrodes (mV)')
        xlabel('time (s)') 
        str = sprintf('Stimulation pair %s for %s', stimchans{stimp,1} , dataBase.NmbrofStims);
        title(str)
        counter = 0;
        nmbr_of_ones = TotOnesStim(stimp);
        name = cell(1,nmbr_of_ones);
    
        % Je wilt weten welke elektroden bij dit stimulatiepaar horen
        % Alleen van die elektroden neem je dan de ccep_plot        
        % Ook alleen die elektroden plot je. 
        
        % Ergens klopt het nog niet met het tellen van het aantal ones, want een aantal elektroden hebben nog niet de juiste naam.
        
        for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)            % for all electrodes

            for i = 1:size(LocaOnes)     
                if ismember({Stimpnm}, LocaOnes{i,1}) && ismember(dataBase.ch{elec}, LocaOnes(i,2)) 
                    counter = counter+1;
                    name(:,counter) =  [LocaOnes(i,2)];                 % names of electrodes with a 1 on stimulation pair stimp.
                end
            end
        end
        
        
        for nmr_plot = 1:nmbr_of_ones                            % total number of ones for stimulation pair stimp  
            Electodes = [dataBase.ch(:)];
        
            loc_electr = all(ismember(Electodes, name{nmr_plot}),2 ) ; 
            loc_elect(:,nmr_plot) = find(loc_electr) ;              % Number of the location of the electrode with a 1 on stimpair stimp.
        
            ccep_plot = squeeze(dataBase.cc_epoch_sorted_avg(loc_elect(nmr_plot),stimp,:));
            ccep_plot(tt>-0.010 & tt<0.010) = NaN;
            plot(tt, nmr_plot*500+ccep_plot, 'k', 'LineWidth', 2);             % plot(tt, loc_elect*500+ccep_plot, 'k', 'LineWidth', 3);      
            hold on
           
        end
       
        
     set(gca,'YTick',500*(1:counter),'YTickLabel',name(:)) ;
   

    % Save the figures
        if dataBase.save_fig==1
            % create folder to save figures
            if ~ exist(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.NmbrofStims),'dir')

                mkdir(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.NmbrofStims));
            end

            % filename
            figureName = fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.NmbrofStims,...
                [dataBase.sub_label '_' dataBase.NmbrofStims '_stimp' stimchans{stimp}]);
            set(gcf,'PaperPositionMode','auto');
            print('-dpng','-r300',figureName);
        else
            pause
        end
    close all
   
               
    end
end

          
