function  scatter_networkPar(dataBase)


for i = 1:size(dataBase,2)
      
    figure('Position',[302,229,938,827])
    subplot(3,1,1)
    scatter(dataBase(i).agreement_parameter.indegreeN_2  ,dataBase(i).agreement_parameter.indegreeN_10  )
    xlabel('2 stims')
    ylabel('10 stims')
    title(sprintf('indegree normalised, %s',dataBase(i).sub_label))
    h = refline;
    h.LineWidth = 2;

    subplot(3,1,2)
    scatter(dataBase(i).agreement_parameter.outdegreeN_2  ,dataBase(i).agreement_parameter.outdegreeN_10  )
    xlabel('2 stims')
    ylabel('10 stims')
    title(sprintf('outdegree normalised, %s',dataBase(i).sub_label))
    h = refline;
    h.LineWidth = 2;

    subplot(3,1,3)
    scatter(dataBase(i).agreement_parameter.BCN_2  ,dataBase(i).agreement_parameter.BCN_10  )
    xlabel('2 stims')
    ylabel('10 stims')
    title(sprintf('BC normalised, %s',dataBase(i).sub_label))
    h = refline;
    h.LineWidth = 2;
    
end

end
