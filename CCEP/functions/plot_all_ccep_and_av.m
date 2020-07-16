function plot_all_ccep_and_av(dataBase,dataBase2, myDataPath, LocOnes, stimchans, dif_mat, TotOnesStim)

tic;
tt = dataBase.tt;    
LocaOnes = LocOnes{:,:}; 

% ER in 10 stims = 1, ER in 2 stims = -1 (dif_mat = 10stim - 2stim --> ER in 10stim (1) - non-ER in 2 stim(0) = 1)
[ER_in10(:,2), ER_in10(:,1)] = find(dif_mat == 1);            % [stimpairs electrodes]

for i = 1:size(ER_in10(:,2))                            % For the number of ones detected
    ER_in10st(i,1) = stimchans(ER_in10(i,1))';
    ER_in10st(i,2) = dataBase.ch(ER_in10(i,2));    
end


set(groot,'defaultFigureVisible','off') % 'on' to turn figures showing on, 'off' to not show the figures.  

[indivstimp,~,stimprow] = unique(sort(dataBase.cc_stimsets,2),'rows');


for stimp = 1:length(indivstimp)                            %1:max(indivstimp)
        
        Stimpnm = stimchans{stimp};     
        stimpnm = find(stimprow == stimp);

       for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)    % for every electrode
            elecnm = dataBase.ch{elec};
            
            ccep_plot = zeros(length(tt),1);
            ccep_plot2 = zeros(length(tt),1); 
            counter_10 = 0;
            nmbr_of_ones = TotOnesStim(stimp);
            NameER = cell(1,nmbr_of_ones);
            
            for i = 1:size(LocaOnes)     
                if ismember({Stimpnm}, LocaOnes{i,1}) && ismember(dataBase.ch{elec}, LocaOnes(i,2)) 

                    % Plot all 10 stimulations
                    test = squeeze(dataBase.cc_epoch_sorted(elec,:,stimpnm,:));             % stimpnm is [1;2]!
                    test2 = reshape(test, size(test,1)*size(test,2), size(test,3));
                    test2(:,tt>-0.01 & tt<0.02) = NaN;
                    
                    % Plot the average of the 10 or 2 stimulations per electrode
                    % cc_epoch_sorted_avg(:,ll,:) =  squeeze(nanmean(squeeze(nanmean(cc_epoch_sorted(:,avg_stim,IC==ll,:),2)),2));

                    ccep_plot = squeeze(dataBase.cc_epoch_sorted_avg(elec,stimp,:));
                    ccep_plot2 = squeeze(dataBase2.cc_epoch_sorted_avg(elec,stimp,:));       % stimulations of the 2 stims
                    ccep_plot(tt>-0.010 & tt<0.010) = NaN;
                    ccep_plot2(tt>-0.010 & tt<0.010) = NaN;
                                        
                
                    figure('Position',[1100 0 700 700])
                    subplot(3,1,1)
                    plot(tt,test2'+[0:500:size(test2,1)*500-1]);                    
                    str = sprintf('ER in 2 stims. Stimulation pair %s for electrode %s', Stimpnm, elecnm);
                    title(str)
                    set(gca,'YTick',2500,'YTickLabel',elecnm) ;                  
                    xlim([-.2 1.5])
                    ylabel(sprintf('10 stimuli for %s', Stimpnm) )
                    xlabel('time (s)') 
                                  
                    subplot(3,1,2)
                    plot(tt,ccep_plot);    % size(test2??)                 
                    str = sprintf('Average of 10 stims for %s and %s', Stimpnm, elecnm);
                    title(str)
                    xlim([-.2 1.5])
                    ylabel('Average (mV)')
                    xlabel('time (s)') 
                    set(gca,'YTick',500,'YTickLabel',elecnm) ;
                                                    
                    
                    subplot(3,1,3)
                    plot(tt,ccep_plot2);    % size(test2??)                
                    str = sprintf('Average of 2 stims for %s and %s', Stimpnm, elecnm);
                    title(str)
                    xlim([-.2 1.5])
                    ylabel('Average (mV)')
                    xlabel('time (s)') 
                    set(gca,'YTick',500,'YTickLabel',elecnm) ;
                    
                    % Standard is an ER evoked in 2 stims and not in 10
                    % stims (random choise)
                    str_main = sprintf('ER evoked by 2 stims in %s',elecnm);        
                    sgtitle(str_main) 
                
                 for j = 1:length(ER_in10st)
                    if ismember({Stimpnm}, ER_in10st{j,1}) && ismember(dataBase.ch{elec}, ER_in10st(j,2))     
                        str_main = sprintf('ER evoked by 10 stims in %s',elecnm);        
                        sgtitle(str_main)
                                    
                    end
                 end  

                 
                  

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
     end
            
            
       end
        close all          
   end
    toc
end

          
