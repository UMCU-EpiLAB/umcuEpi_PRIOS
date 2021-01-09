function  scatter_ranking(dataBase,myDataPath)

mode = {'ERs evoked per stimulation pair','Indegree','Outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)

   fig = figure('Position',[302,17,1110,1039]);

    for i = 1:size(dataBase,2)
       if strcmp(mode{J},'ERs evoked per stimulation pair')
            clin = dataBase(i).rank.unsort_SPES_clin(:,4);
            prop =  dataBase(i).rank.unsort_SPES_prop(:,4);
            idx_nan = isnan(dataBase(i).rank.unsort_SPES_clin(:,3)) | isnan(dataBase(i).rank.unsort_SPES_prop(:,3));        % Find NaN values to avoid plotting these           
            pval = dataBase(i).statistics.p_stimp  ;
            rho = dataBase(i).statistics.rho_stimp;

        elseif strcmp(mode{J},'Indegree')
            clin = dataBase(i).rank.unsort_indegreeSPES_clin(:,3);
            prop = dataBase(i).rank.unsort_indegreeSPES_prop(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_indegreeSPES_clin(:,2)) | isnan(dataBase(i).rank.unsort_indegreeSPES_prop(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_indegree;
            rho = dataBase(i).statistics.rho_indegree;

        elseif strcmp(mode{J},'Outdegree')
            clin = dataBase(i).rank.unsort_outdegreeSPES_clin(:,3);
            prop = dataBase(i).rank.unsort_outdegreeSPES_prop(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_outdegreeSPES_clin(:,2)) | isnan(dataBase(i).rank.unsort_outdegreeSPES_prop(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_outdegree;
            rho = dataBase(i).statistics.rho_outdegree;

        elseif strcmp(mode{J},'Betweenness Centrality')
            clin = dataBase(i).rank.unsort_BCSPES_clin(:,3);
            prop = dataBase(i).rank.unsort_BCSPES_prop(:,3);
            idx_nan = isnan(dataBase(i).rank.unsort_BCSPES_clin(:,2)) | isnan(dataBase(i).rank.unsort_BCSPES_prop(:,2));        % Find NaN values to avoid plotting these
            pval = dataBase(i).statistics.p_BC;
            rho = dataBase(i).statistics.rho_BC;
       end
       
       
        subplot(size(dataBase,2),1,i);
        scatter(clin(~idx_nan)  , prop(~idx_nan) ,'*' )
        % Change fontsize
        ax = gca;
        ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval))
    
        xmin = min(clin(~idx_nan));
        xmax = max(clin(~idx_nan)); 
           
        if pval < 0.05
            [P,S] = polyfit(clin(~idx_nan),prop(~idx_nan),1);
            
            [y_fit, delta] = polyval(P,clin(~idx_nan),S);
            
            plot(clin(~idx_nan),prop(~idx_nan),'*')                        % This is equal to scatter
            hold on
            
            % Plot polyfit throught data points
            plot(clin(~idx_nan),y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
            hold on
            % Change fontsize
            ax = gca;
            ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;                      

            % Plot conficence interval as a line
%             plot(clin, y_fit-2*delta, 'm--', clin, y_fit+2*delta,'--','color',[0.6,0.1,0.2,0.8])

            if pval < 0.01
                    title(sprintf('%s, p = <0.01, r_s = %1.3f', dataBase(i).sub_label, rho));
             elseif pval<0.05
                     title(sprintf('%s, p = <0.05, r_s = %1.3f', dataBase(i).sub_label, rho));
            end
            
            % Plot confidence interval as a patch
            Filled_CI = patch([min(clin),max(clin),max(clin),min(clin)], [min(y_fit-2*delta),max(y_fit-2*delta),max(y_fit+2*delta), min(y_fit+2*delta)],[0.1,0.2,0.2]);
            alpha(0.06)                % set patches transparency to 0.
            Filled_CI.EdgeAlpha = 0;
           
        end
        xlim([xmin xmax])
       
   % Create main title without decreasing the subplots sizes
    annotation(gcf,'textbox',[0.32 0.95 0.35 0.043],'VerticalAlignment','middle','String',sprintf('%s',mode{J}),'HorizontalAlignment','center','FontWeight','bold','FontSize',15,'FitBoxToText','off','EdgeColor','none');

   %Create one main x- and y-label
    han=axes(fig,'visible','off'); 
    han.YLabel.Visible='on';
    han.YLabel.Position = [-0.0464    0.5000         0];
    ylabel(han,{'Rank SPES-prop'}, 'FontSize',17,'fontweight','bold')
    
    han.XLabel.Visible='on';
    xlabel(han,{'Rank SPES-clin'}, 'FontSize',17,'fontweight','bold')


    end
    
    % Save figure
    outlabel=sprintf('All_pat_scatter_%s.jpg',mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/Scatter of Ranking/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')    

end   

end
         