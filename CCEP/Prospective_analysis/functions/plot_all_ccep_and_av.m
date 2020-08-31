function plot_all_ccep_and_av(dataBase, dataBase2, myDataPath, LocOnes, agreement)

tic;
dif_mat = agreement.dif_mat;
% TotOnesStim = agreement.TotOnesStim;
tt = dataBase.tt;
LocaOnes = LocOnes{:,:};

% ER in 10 stims = 1, ER in 2 stims = -1 (dif_mat = 10stim - 2stim --> ER in 10stim (1) - non-ER in 2 stim(0) = 1)
[ER_in10(:,2), ER_in10(:,1)] = find(dif_mat == 1);            % [stimpairs electrodes]

ER_in10st = cell(size(ER_in10));
for i = 1:size(ER_in10(:,2))                            % For the number of ones detected
    ER_in10st(i,1) = dataBase.stimpnames_avg(ER_in10(i,1))';
    ER_in10st(i,2) = dataBase.ch(ER_in10(i,2));
end

set(groot,'defaultFigureVisible','on') % 'on' to turn figures showing on, 'off' to not show the figures.

% [indivstimp,~,stimprow] = unique(sort(dataBase.cc_stimsets_all,2),'rows');

for stimp = 1:size(dataBase.stimsets_avg,1)                           %1:max(indivstimp)
    
    Stimpnm = dataBase.stimpnames_avg{stimp};
    %     stimpnm = find(stimprow == stimp);
    
    for elec = 1:size(dataBase.ch)    % for every electrode
        elecnm = dataBase.ch{elec};
        
        %         ccep_plot = zeros(length(tt),1);
        %         ccep_plot2 = zeros(length(tt),1);
        %         counter_10 = 0;
%                 nmbr_of_ones = TotOnesStim(stimp);
        %         NameER = cell(1,nmbr_of_ones);
        
        %         for i = 1:size(LocaOnes)
        if any(strcmp(LocaOnes(:,1),Stimpnm) & strcmp(LocaOnes(:,2),elecnm))        % Check if the stimpair has evoked an CCEP in each electrode
            
            % Plot all 10 stimulations
            ccep_plot = squeeze(dataBase.epoch_sorted_select_avg(elec,stimp,:,:));             % stimpnm is [1;2]!
            ccep_plot(:,tt>-0.01 & tt<0.02) = NaN;
            ccep_plot2 = squeeze(dataBase2.epoch_sorted_select_avg(elec,stimp,:,:));             % stimpnm is [1;2]!
            ccep_plot2(:,tt>-0.01 & tt<0.02) = NaN;
            
            % Plot the average of the 10 or 2 stimulations per electrode
            % cc_epoch_sorted_avg(:,ll,:) =  squeeze(nanmean(squeeze(nanmean(cc_epoch_sorted(:,avg_stim,IC==ll,:),2)),2));
            
            ccep_plot_avg = squeeze(dataBase.epoch_sorted_avg(elec,stimp,:));
            ccep_plot2_avg = squeeze(dataBase2.epoch_sorted_avg(elec,stimp,:));       % stimulations of the 2 stims
            ccep_plot_avg(tt>-0.010 & tt<0.010) = NaN;
            ccep_plot2_avg(tt>-0.010 & tt<0.010) = NaN;
            
            figure('Position',[1100 0 700 700])
            subplot(3,1,1)
            plot(tt,ccep_plot'+(0:500:size(ccep_plot,1)*500-1));
            str = sprintf('Stimulation pair %s for electrode %s', Stimpnm, elecnm);
            title(str)
            set(gca,'YTick',2500,'YTickLabel',elecnm) ;
            xlim([-.2 1.5])
            ylabel(sprintf('10 stimuli for %s', Stimpnm) )
            xlabel('time (s)')
            
            subplot(3,1,2)
            plot(tt,ccep_plot','Color',[0.7 0.7 0.7]);
            hold on
            plot(tt,ccep_plot_avg,'k','LineWidth',2);
            hold off
            str = sprintf('Average of 10 stims for %s and %s', Stimpnm, elecnm);
            title(str)
            xlim([-.2 1.5])
            ylabel('Average (mV)')
            xlabel('time (s)')
            set(gca,'YTick',500,'YTickLabel',elecnm) ;
            
            subplot(3,1,3)
            plot(tt,ccep_plot2','Color',[0.7 0.7 0.7]);
            hold on
            plot(tt,ccep_plot2_avg,'k','LineWidth',2);
            hold off
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
            
            if any(strcmp(ER_in10st(:,1),Stimpnm) & strcmp(ER_in10st(:,2),elecnm))
                str_main = sprintf('ER evoked by 10 stims in %s',elecnm);
                sgtitle(str_main)
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
        close all
    end
end
toc
end


