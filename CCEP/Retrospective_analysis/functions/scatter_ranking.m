function  scatter_ranking(dataBase,myDataPath)

mode = {'ERs evoked per stimulation pair','Indegree','Outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)

   figure('Position',[302,17,1110,1039])


    for i = 1:size(dataBase,2)
       if strcmp(mode{J},'ERs evoked per stimulation pair')
            par10 = dataBase(i).rank.unsort_stims10(:,4);
            par2 =  dataBase(i).rank.unsort_stims2(:,4);
            idx_nan = isnan(dataBase(i).rank.unsort_stims10(:,3)) | isnan(dataBase(i).rank.unsort_stims2(:,3));        % Find NaN values to avoid plotting these           
            pval = dataBase(i).statistics.p_stimp  ;
            rho = dataBase(i).statistics.rho_stimp;

        elseif strcmp(mode{J},'Indegree')
            par10 = dataBase(i).rank.unsort_indegreestims10(:,3);
            par2 = dataBase(i).rank.unsort_indegreestims2(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_indegreestims10(:,2)) | isnan(dataBase(i).rank.unsort_indegreestims2(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_indegree;
            rho = dataBase(i).statistics.rho_indegree;

        elseif strcmp(mode{J},'Outdegree')
            par10 = dataBase(i).rank.unsort_outdegreestims10(:,3);
            par2 = dataBase(i).rank.unsort_outdegreestims2(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_outdegreestims10(:,2)) | isnan(dataBase(i).rank.unsort_outdegreestims2(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_outdegree;
            rho = dataBase(i).statistics.rho_outdegree;

        elseif strcmp(mode{J},'Betweenness Centrality')
            par10 = dataBase(i).rank.unsort_BCstims10(:,3);
            par2 = dataBase(i).rank.unsort_BCstims2(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_BCstims10(:,2)) | isnan(dataBase(i).rank.unsort_BCstims2(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_BC;
            rho = dataBase(i).statistics.rho_BC;
       end
       
       

        subplot(size(dataBase,2),1,i);
        scatter(par10(~idx_nan)  , par2(~idx_nan) ,'*' )
        ylabel({'Place on ranking';'2 stims'}, 'FontSize',11)
        xlabel("Place on ranking 10 stims",'FontSize',11)
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval))
%         legend(sprintf('%s',mode{J}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',11)
%         ax = gca;
%         ax.YAxis.FontSize = 10;
%         ax.XAxis.FontSize = 10;     
        xmin = min(par10);
        xmax = max(par10(~idx_nan));
                      
           
        if pval < 0.05
            [P,S] = polyfit(par10(~idx_nan),par2(~idx_nan),1);
            
            [y_fit, delta] = polyval(P,par10(~idx_nan),S);
            
            plot(par10(~idx_nan),par2(~idx_nan),'*')                        % This is equal to scatter
            hold on
            
            % Plot polyfit throught data points
            plot(par10(~idx_nan),y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
            hold on
            ylabel({'Place on ranking';'2 stims'}, 'FontSize',11)
            xlabel("Place on ranking 10 stims",'FontSize',11)
             
            if pval < 0.01
                title(sprintf('%s, p = <0.01, r_s = %1.3f', dataBase(i).sub_label, rho));
            elseif pval<0.05
                 title(sprintf('%s, p = <0.05, r_s = %1.3f', dataBase(i).sub_label, rho));
            end
            % Plot conficence interval as a line
%             plot(par10, y_fit-2*delta, 'm--', par10, y_fit+2*delta,'--','color',[0.6,0.1,0.2,0.8])
            
            % Plot confidence interval as a patch
            Filled_CI = patch([min(par10),max(par10),max(par10),min(par10)], [min(y_fit-2*delta),max(y_fit-2*delta),max(y_fit+2*delta), min(y_fit+2*delta)],[0.1,0.2,0.2]);
            alpha(0.06)                % set patches transparency to 0.
            Filled_CI.EdgeAlpha = 0;
%             title(sprintf('%s, p = %1.3f', dataBase(i).sub_label,pval));
           
        end
        xlim([xmin xmax])
       
       % Create main title without decreasing the subplots sizes
        annotation(gcf,'textbox',[0.32 0.95 0.35 0.043],'VerticalAlignment','middle','String',sprintf('%s',mode{J}),'HorizontalAlignment','center','FontWeight','bold','FontSize',15,'FitBoxToText','off','EdgeColor','none');
        
    end
    
   % Save figure
    outlabel=sprintf('All_pat_scatter_%s.jpg',mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/Scatter van ranking/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')    

end   

end
         