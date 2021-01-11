function  scatter_networkPar(dataBase, myDataPath)

mode = {'ERs per stimulation pair','Indegree','Outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)
    
    fig = figure('Position',[302,17,1224,1039]);
    
    for i = 1:size(dataBase,2)
        
        if strcmp(mode{J},'ERs per stimulation pair')
            par10 = dataBase(i).agreement_parameter.ERs_stimp10;
            par2 = dataBase(i).agreement_parameter.ERs_stimp2;
            pval = dataBase(i).statistics.p_stimp;
            rho = dataBase(i).statistics.rho_stimp;
            
        elseif strcmp(mode{J},'Indegree')
            par10 = dataBase(i).agreement_parameter.indegreeN_10;
            par2 = dataBase(i).agreement_parameter.indegreeN_2;
            pval = dataBase(i).statistics.p_indegree;
            rho = dataBase(i).statistics.rho_indegree;
            
        elseif strcmp(mode{J},'Outdegree')
            par10 = dataBase(i).agreement_parameter.outdegreeN_10;
            par2 = dataBase(i).agreement_parameter.outdegreeN_2;
            pval = dataBase(i).statistics.p_outdegree;
            rho = dataBase(i).statistics.rho_outdegree;
            
        elseif strcmp(mode{J},'Betweenness Centrality')
            par10 = dataBase(i).agreement_parameter.BCN_10;
            par2 = dataBase(i).agreement_parameter.BCN_2;
            pval = dataBase(i).statistics.p_BC;
            rho = dataBase(i).statistics.rho_BC;
        end
        
        % Find NaN values to avoid plotting these
        idx_nan = isnan(par10) | isnan(par2);   
        
        subplot(size(dataBase,2),1,i)
        scatter(par10  , par2 ,'*' )
        
        % Change fontsize
        ax = gca;
        ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval),'FontSize',12)
        xmin = min(par10(~idx_nan));
        xmax = max(par10(~idx_nan));
        
        
        if pval < 0.05
            [P,S] = polyfit(par10(~idx_nan),par2(~idx_nan),1);
            [y_fit, ~] = polyval(P,par10(~idx_nan),S);
            
            plot(par10(~idx_nan),par2(~idx_nan),'*')                        % This is equal to scatter
            hold on
            
            % Plot polyfit throught data points
            plot(par10(~idx_nan),y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
            hold on
            % Change fontsize
            ax = gca;
            ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;
            
            % Plot conficence interval as a line
%             plot(par10, y_fit-2*delta, 'm--', par10, y_fit+2*delta,'--','color',[0.6,0.1,0.2,0.8])
            
%             % Plot confidence interval as a patch
%             Filled_CI = patch([min(par10),max(par10),max(par10),min(par10)], [min(y_fit-2*delta),max(y_fit-2*delta),max(y_fit+2*delta), min(y_fit+2*delta)],[0.1,0.2,0.2]);
%             alpha(0.06)                % set patches transparency to 0.
%             Filled_CI.EdgeAlpha = 0;
%             title(sprintf('%s, p = %1.3f', dataBase(i).sub_label,pval));
%             
            if pval < 0.01
                title(sprintf('%s, p = <0.01, r_s = %1.3f', dataBase(i).sub_label, rho),'FontSize',12)
            elseif pval<0.05
                title(sprintf('%s, p = <0.05, r_s = %1.3f', dataBase(i).sub_label, rho),'FontSize',12)
            end
        
             
%             P = polyfit(par10(~idx_nan),par2(~idx_nan),1);
%             X = xmin : 0.1*xmax : xmax+0.2*xmax;
%             Y = P(1)*X + P(2);
%             
%             hold on
%             h=plot(X,Y);
%             hold off
%             h.LineWidth = 2;
%             title(sprintf('%s, p = %1.3f', dataBase(i).sub_label,pval))
% %             legend(sprintf('%s',mode{J}), sprintf('r_s = %1.3f',rho  ),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)
%             
            
        end
        xlim([xmin xmax])
        
        % Create main title without decreasing the subplots sizes
        annotation(gcf,'textbox',[0.32 0.95 0.35 0.043],'VerticalAlignment','middle','String',sprintf('%s',mode{J}),'HorizontalAlignment','center','FontWeight','bold','FontSize',15,'FitBoxToText','off','EdgeColor','none');
       
        %Create one main x- and y-label
        han=axes(fig,'visible','off'); 
        han.YLabel.Visible='on';
        han.YLabel.Position = [-0.0864    0.5000         0];
        ylabel(han,{'Value 2 stimuli setting'}, 'FontSize',17,'fontweight','bold')
    
        han.XLabel.Visible='on';
        xlabel(han,{'Value 10 stimuli setting'}, 'FontSize',17,'fontweight','bold')

    end
    
    % Save figure
    outlabel=sprintf('All_pat_scatter_%s.jpg',mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')
    
end

end


