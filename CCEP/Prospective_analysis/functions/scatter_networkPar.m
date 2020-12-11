function  scatter_networkPar(dataBase,myDataPath)

mode = {'ERs per stimulation pair','Indegree','Outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)

   figure('Position',[302,17,1110,1039])


    for i = 1:size(dataBase,2)
       if strcmp(mode{J},'ERs per stimulation pair')
            clin = dataBase(i).agreement_parameter.ERs_stimpClin;
            prop = dataBase(i).agreement_parameter.ERs_stimpProp;
            pval = dataBase(i).statistics.p_stimp;
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
        
        subplot(size(dataBase,2),1,i);
        scatter(clin  , prop ,'*' )
        ylabel('value SPES-prop', 'FontSize',9)
        xlabel("Value for SPES-Clin"+newline+"    ",'FontSize',9)
%         str_main = sprintf('%s', mode{J});
%         sgtitle(str_main)
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval))
        legend(sprintf('%s',mode{J}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',11)
        ax = gca;
        ax.YAxis.FontSize = 10;
        ax.XAxis.FontSize = 10;     
        xmin = 0;
        xmax = round(max(clin)+0.1*max(clin),2);
        
        if pval < 0.05
            idx_nan = isnan(clin) | isnan(prop);
            P = polyfit(clin(~idx_nan),prop(~idx_nan),1);
            X = xmin:0.1*xmax:xmax+0.2*xmax;
            Y = P(1)*X + P(2);
            
            hold on
            h=plot(X,Y,'LineWidth',2);
            hold off
            h.Color(4) = 0.7;
            uistack(h,'bottom')         %you can also do uistack(h,'top')
            title(sprintf('%s, p = <0.05', dataBase(i).sub_label))
            legend(sprintf('%s',mode{J}), sprintf('r_s = %1.3f',rho  ),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',11)
            if pval < 0.01
                title(sprintf('%s, p = <0.01', dataBase(i).sub_label))
            end
            
        end
        xlim([xmin xmax])
        
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
         