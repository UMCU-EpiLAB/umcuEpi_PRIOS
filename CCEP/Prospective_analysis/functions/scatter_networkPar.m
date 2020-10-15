function  scatter_networkPar(dataBase,myDataPath)


for i = 1:size(dataBase,2)
      
    figure('Position',[302,17,938,1039])
    subplot(4,1,1)
    scatter(dataBase(i).agreement_parameter.ERs_stimpProp  ,dataBase(i).agreement_parameter.ERs_stimpClin  )
    xlabel('SPES-prop')
    ylabel('SPES-clin')
    title(sprintf('ERs per stimulation pair, %s, p =  %1.3f',dataBase(i).sub_label, dataBase(i).statistics.p_stimp))
    legend('ERs per stimulation', sprintf('rho = %1.3f',dataBase(i).statistics.rho_stimp  ));
    
    if dataBase(i).statistics.p_stimp < 0.05
        h = refline;
        h.LineWidth = 2;
    end
        
    subplot(4,1,2)
    scatter(dataBase(i).agreement_parameter.indegreeN_Prop  ,dataBase(i).agreement_parameter.indegreeN_Clin  )
    xlabel('SPES-prop')
    ylabel('SPES-clin')
    title(sprintf('indegree normalised, %s, p = %1.3f',dataBase(i).sub_label, dataBase(i).statistics.p_indegree))
    h = refline;
    h.LineWidth = 2;
    legend('Intdegree', sprintf('rho = %1.3f',dataBase(i).statistics.rho_indegree  ));
     if dataBase(i).statistics.p_indegree < 0.05
        h = refline;
        h.LineWidth = 2;
    end

    subplot(4,1,3)
    scatter(dataBase(i).agreement_parameter.outdegreeN_Prop  ,dataBase(i).agreement_parameter.outdegreeN_Clin  )
    xlabel('SPES-prop')
    ylabel('SPES-clin')
    title(sprintf('outdegree normalised, %s, p =  %1.3f',dataBase(i).sub_label,  dataBase(i).statistics.p_outdegree))
    h = refline;
    h.LineWidth = 2;
    legend('Outdegree', sprintf('rho = %1.3f',dataBase(i).statistics.rho_outdegree  ));
     if dataBase(i).statistics.p_outdegree < 0.05
        h = refline;
        h.LineWidth = 2;
    end


    subplot(4,1,4)
    scatter(dataBase(i).agreement_parameter.BCN_Prop  ,dataBase(i).agreement_parameter.BCN_Clin  )
    xlabel('SPES-prop')
    ylabel('SPES-clin')
    title(sprintf('BC normalised, %s, p = %1.3f',dataBase(i).sub_label, dataBase(i).statistics.p_BC))
    h = refline;
    h.LineWidth = 2;
    legend('BC', sprintf('rho = %1.3f',dataBase(i).statistics.rho_BC  ));
    if dataBase(i).statistics.p_BC < 0.05
        h = refline;
        h.LineWidth = 2;
    end
    
    
    % Save figure
    outlabel=sprintf('sub-%s_scatterNetwPar.jpg',dataBase(i).sub_label);
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')


end

end
