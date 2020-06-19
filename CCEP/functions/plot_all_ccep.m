function plot_all_ccep(dataBase, myDataPath, stimchans)
tic;
ccep_plot = zeros(10,8192);
tt = dataBase.tt;    
   for stimp = 1:2:size(dataBase.cc_epoch_sorted,3)                          % for all stimulation pairs
       % uneven stimpnumbers are positve direction
       % even stimpnumbers are negative direction
       
      
       
       for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)                     % for every electrode
           
           for j = 1:2:2*(size(dataBase.cc_epoch_sorted,2))                  % for every stimulation
                 ccep_plot(j:j+1,:,:) = squeeze(dataBase.cc_epoch_sorted(elec,round(j/2),[stimp:(stimp+1)],:));     % all stimuli of the positive and negative direction
                 ccep_plot(tt>-0.010 & tt<0.010) = NaN;                      % Make signal around artefact NaN
           end
       
        figure('Position',[200 0 700 700])
        plot(tt, 4500 +ccep_plot(1,:,:), 'LineWidth',3)                 % this is stimulation 1 of positive direction
        hold on
        plot(tt, 4000 +ccep_plot(3,:,:))
        plot(tt, 3500 +ccep_plot(5,:,:))
        plot(tt, 3000 +ccep_plot(7,:,:))
        plot(tt, 2500 +ccep_plot(9,:,:))

        plot(tt, 2000 + ccep_plot(2,:,:), 'LineWidth',3)                 % this are stimulation 1 of negative direction
        plot(tt, 1500 +ccep_plot(4,:,:))
        plot(tt, 1000 +ccep_plot(6,:,:))
        plot(tt, 500 +ccep_plot(8,:,:))
        plot(tt, ccep_plot(10,:,:))
        stimuli_nm = {'neg5' 'neg4' 'neg3' 'neg2' 'neg1' 'pos5' 'pos4' 'pos3' 'pos2' 'pos1'};
        legend('pos1','pos2','pos3','pos4','pos5','neg1','neg2','neg3','neg4','neg5')
               
        set(gca,'YTick',500*(0:size(ccep_plot)-1), 'YTickLabel',stimuli_nm) 
        Stimpnm = dataBase.stimpnames{stimp};
        elecnm = dataBase.ch{elec};
        xlim([-.2 1.5])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 
        str = sprintf('Stimulation pair %s for electrode %s', Stimpnm, elecnm);
        title(str)
        hold off

        % Save the figures
            if dataBase.save_fig==1
                % create folder to save figures
                if ~ exist(fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm),'dir')

                    mkdir(fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm));
                end

                % filename
                figureName = fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm,...
                    [dataBase.sub_label '_stimp_' Stimpnm '_elec_' elecnm ]);
                set(gcf,'PaperPositionMode','auto');
                print('-dpng','-r300',figureName);
            else
                pause
            end
       end
        close all          
   end
    toc
end

          
