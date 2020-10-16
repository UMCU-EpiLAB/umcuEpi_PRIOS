function  scatter_networkPar(dataBase, myDataPath)

mode = {'ERs per stimulation pair','indegree','outdegree','Betweenness Centrality'};

for J = 1:size(mode,2)
    
    figure('Position',[302,17,938,1039])

    for i = 1:size(dataBase,2)
        
    if strcmp(mode{J},'ERs per stimulation pair')
        par10 = dataBase(i).agreement_parameter.ERs_stimp10;
        par2 = dataBase(i).agreement_parameter.ERs_stimp2;
        pval = dataBase(i).statistics.p_stimp;
        rho = dataBase(i).statistics.rho_stimp;
        
    elseif strcmp(mode{J},'indegree')
        par10 = dataBase(i).agreement_parameter.indegreeN_10;
        par2 = dataBase(i).agreement_parameter.indegreeN_2;
        pval = dataBase(i).statistics.p_indegree;
        rho = dataBase(i).statistics.rho_indegree;
        
    elseif strcmp(mode{J},'outdegree')
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
    
    
        subplot(size(dataBase,2),1,i)
        scatter(par2  , par10  )
        ylabel('10 stimuli situation')
        xlabel("2 stimuli situation"+newline+"   ")
        str_main = sprintf('%s', mode{J});
        sgtitle(str_main)
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval))
        legend(sprintf('%s',mode{J}))
      
        if pval < 0.05
            h = refline;
            h.LineWidth = 2;
            legend(sprintf('%s',mode{J}), sprintf('rho = %1.3f',rho  ));
        end
        hold on

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


