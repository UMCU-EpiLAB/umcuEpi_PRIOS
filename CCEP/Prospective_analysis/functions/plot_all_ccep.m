function plot_all_ccep(ccep_data, myDataPath)
%%% KLOPT NOG NIET!!! HIJ LAAT NU BIJ EL01-EL02 NIET DE FOUT ZIEN ZOALS
%%% NORMAAL BIJ EL01 DAN TE ZIEN MOET ZIJN. 


tt = ccep_data.tt;    
set(groot,'defaultFigureVisible','off') % 'on' to turn figures showing on, 'off' to not show the figures.  

% When the status of the channel is bad, not visualise the response
bad_sig(:,1) = find(contains(ccep_data.tb_channels.status, {'bad'}));  
ccep_data.cc_epoch_sorted_avg(bad_sig,:,:) = NaN;


    for stimp = 1:length(ccep_data.cc_stimsets_avg)            % For each stimulation pair                              
        Stimpnm = ccep_data.stimpnames_avg{stimp};     
        
        % Remove electrodes which are stimulated
        stim_elec = ccep_data.cc_stimsets_avg(stimp,:);
        plot_ccep = ccep_data.cc_epoch_sorted_avg(:,:,:);
        plot_ccep(stim_elec,:,:) = NaN;

        
        ccep = squeeze(plot_ccep(:,stimp,:));
        ccep(:,tt>-0.01 & tt<0.02) = NaN;

        figure('Position',[896,21,1004,1041])
        plot(tt, ccep' + [0:1000:size(ccep,1)*1000-1])

        str = sprintf('%s', Stimpnm);
        title(str)
        set(gca,'YTick',[0:1000:size(ccep,1)*1000-1],'YTickLabel',ccep_data.ch) ;                  
        xlim([-.2 1.5])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 

       % Save the figures
            if ccep_data.save_fig=='y'
                % create folder to save figures
                if ~ exist(fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_data.sub_label,ccep_data.task_name),'dir')

                    mkdir(fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_data.sub_label,ccep_data.task_name));
                end

                % filename
                figureName = fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_data.sub_label,ccep_data.task_name,...
                    [ccep_data.sub_label '_stimp_' Stimpnm ]);
                set(gcf,'PaperPositionMode','auto');
                print('-dpng','-r300',figureName);
            else
                pause
            end
            
    end
close all
end
            

          
