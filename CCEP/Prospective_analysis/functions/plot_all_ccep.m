function plot_all_ccep(ccep_clin, myDataPath)
%%% KLOPT NOG NIET!!! HIJ LAAT NU BIJ EL01-EL02 NIET DE FOUT ZIEN ZOALS
%%% NORMAAL BIJ EL01 DAN TE ZIEN MOET ZIJN. 


tt = ccep_clin.tt;    
set(groot,'defaultFigureVisible','off') % 'on' to turn figures showing on, 'off' to not show the figures.  

for stimp = 1:length(ccep_clin.cc_stimsets_avg)            % For each stimulation pair                           
           
         Stimpnm = ccep_clin.stimpnames_avg{stimp};     

       for elec = 1:size(ccep_clin.ch,1)                     % for every electrode
             elecnm = ccep_clin.ch{elec};

                    ccep = squeeze(ccep_clin.cc_epoch_sorted(elec,:,stimp,:));
                    ccep(:,tt>-0.01 & tt<0.02) = NaN;
                
                    figure('Position',[1100 0 700 700])
                    plot(tt, ccep'+ [0:1000:size(ccep,1)*1000-1])
                    
                    str = sprintf('%s on %s', Stimpnm, elecnm);
                    title(str)
                    set(gca,'YTick',2500,'YTickLabel',elecnm) ;                  
                    xlim([-.2 1.5])
                    ylabel('All stimuli of this stimulation pair' )
                    xlabel('time (s)') 

        % Save the figures
            if ccep_clin.save_fig=='y'
                % create folder to save figures
                if ~ exist(fullfile(myDataPath.CCEPpath,'all_ccep_figures',ccep_clin.sub_label,Stimpnm),'dir')

                    mkdir(fullfile(myDataPath.CCEPpath,'all_ccep_figures',ccep_clin.sub_label,Stimpnm));
                end

                % filename
                figureName = fullfile(myDataPath.CCEPpath,'all_ccep_figures',ccep_clin.sub_label,Stimpnm,...
                    [ccep_clin.sub_label '_stimp_' Stimpnm '_elec_' elecnm ]);
                set(gcf,'PaperPositionMode','auto');
                print('-dpng','-r300',figureName);
            else
                pause
            end
            
       end
        close all
end
            
end

          
