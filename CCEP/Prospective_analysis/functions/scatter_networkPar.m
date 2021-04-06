function  scatter_networkPar(dataBase,myDataPath)
% Visualise the agreement in a scatter plot of the absolute values
% The significance of the absolute values is determined with the Wilcoxon
% signed rank test and this does not result in an Rho value. 
% When statistically significant --> there if a difference between the 2
% measurements.

mode = {'ERs evoked per stimulation pair','Indegree','Outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)

   fig = figure('Position',[302,17,1110,1039]);

    for i = 1:size(dataBase,2)
       if strcmp(mode{J},'ERs evoked per stimulation pair')
            clin = dataBase(i).agreement_parameter.ERs_stimpClin;
            prop = dataBase(i).agreement_parameter.ERs_stimpProp;
            pval = dataBase(i).statistics.p_stimp  ;
            rho = dataBase(i).statistics.rho_stimp;

        elseif strcmp(mode{J},'Indegree')
            clin = dataBase(i).agreement_parameter.indegreeN_Clin;
            prop = dataBase(i).agreement_parameter.indegreeN_Prop;
            pval = dataBase(i).statistics.p_indegree;
            rho = dataBase(i).statistics.rho_indegree;
            
        elseif strcmp(mode{J},'Outdegree')
            clin = dataBase(i).agreement_parameter.outdegreeN_Clin;
            prop = dataBase(i).agreement_parameter.outdegreeN_Prop;
            pval = dataBase(i).statistics.p_outdegree;
            rho = dataBase(i).statistics.rho_outdegree;

        elseif strcmp(mode{J},'Betweenness Centrality')
            clin = dataBase(i).agreement_parameter.BCN_Clin;
            prop = dataBase(i).agreement_parameter.BCN_Prop;
            pval = dataBase(i).statistics.p_BC;
            rho = dataBase(i).statistics.rho_BC;
       end
       
        % Find NaN values to avoid plotting these
        idx_nan = isnan(clin) | isnan(prop);        
 
        subplot(size(dataBase,2),1,i);
        plot(clin(~idx_nan)  , prop(~idx_nan) ,'*' )
        box off

        % Change fontsize
        ax = gca;
        ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;
           
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval),'FontSize',12)
        xmin = min(clin(~idx_nan));
        xmax = max(clin(~idx_nan));
             
                   
        if pval < 0.05
            [P,S] = polyfit(clin(~idx_nan),prop(~idx_nan),1);
            [y_fit, ~] = polyval(P,clin(~idx_nan),S);
            
            plot(clin(~idx_nan),prop(~idx_nan),'*') 
            box off
            hold on
            
            % Plot polyfit throught data points
            plot(clin(~idx_nan),y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
            hold on
            % Change fontsize
            ax = gca;
            ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;

            % Plot conficence interval as a line
%             plot(clin, y_fit-2*delta, 'm--', clin, y_fit+2*delta,'--','color',[0.6,0.1,0.2,0.8])
            
            % Plot confidence interval as a patch
%             Filled_CI = patch([min(clin),max(clin),max(clin),min(clin)], [min(y_fit-2*delta),max(y_fit-2*delta),max(y_fit+2*delta), min(y_fit+2*delta)],[0.1,0.2,0.2]);
%             alpha(0.06)                % set patches transparency to 0.
%             Filled_CI.EdgeAlpha = 0;
            
            if pval < 0.001
                title(sprintf('%s, p = <0.001, r_s = %1.3f', dataBase(i).sub_label, rho),'FontSize',12)
            else
                 title(sprintf('%s, p = %1.3f, r_s = %1.3f', dataBase(i).sub_label, pval, rho),'FontSize',12)
            end   
            
%             legend(sprintf('%s',mode{J}), sprintf('r_s = %1.3f',rho  ),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',11)
%            X = xmin:0.1*xmax:xmax+0.2*xmax;
%             Y = P(1)*X + P(2);
%             
%             hold on
%             h=plot(X,Y,'LineWidth',2);
%             hold off
%             h.Color(4) = 0.7;
%             uistack(h,'bottom')         %you can also do uistack(h,'top')
           
        end
        xlim([xmin xmax])
             
        % Create main title without decreasing the subplots sizes
        annotation(gcf,'textbox',[0.32 0.95 0.35 0.043],'VerticalAlignment','middle','String',sprintf('%s',mode{J}),'HorizontalAlignment','center','FontWeight','bold','FontSize',15,'FitBoxToText','off','EdgeColor','none');
        
    end
    
     %Create one main x- and y-label
     han=axes(fig,'visible','off'); 
     han.YLabel.Visible='on';
     han.YLabel.Position = [-0.0584    0.5000         0];
     ylabel(han,{'Value SPES-prop'}, 'FontSize',17,'fontweight','bold')
    
     han.XLabel.Visible='on';
     xlabel(han,{'Value SPES-clin'}, 'FontSize',17,'fontweight','bold')

        
    % Save figure
    outlabel=sprintf('All_pat_scatter_%s.png',mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'png')    

end   

end
         