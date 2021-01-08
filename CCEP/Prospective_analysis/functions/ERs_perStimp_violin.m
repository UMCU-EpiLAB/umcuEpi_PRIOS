function ERs_perStimp_violin(dataBase,myDataPath)
% Create a violin plot for the number of ERs evoked per stimulation pair. 
% The stimulation pairs are connected via lines. 
new_mat =NaN([54,12]);
    
for subj = 1:size(dataBase,2)
    Clin = dataBase(subj).agreement_parameter.ERs_stimpClin;
    Prop = dataBase(subj).agreement_parameter.ERs_stimpProp;   
    sub_label = dataBase(subj).sub_label;
        
     i = 1;
     clin_colm = 2*subj-1;                      % prealloction of the column number
     prop_colm = 2*subj;                        % prealloction of the column number

        for stimp = 1:size(Clin,1)                          % For each stimpair                                          
            new_mat(i,clin_colm) = Clin(stimp);            % plot the SPES-clin amp in column 1
            new_mat(i,prop_colm) = Prop(stimp);          % plot the SPES-prop amp in column 2
            i = i+1;           
        end               
end

    
    figure('Position',[205,424,1530,638]);
    grouporder = {'PRIOS01**','','  PRIOS02*','','  PRIOS03','','  PRIOS04**','','  PRIOS05**','','  PRIOS06**',''};

    violins = violinplot([new_mat],grouporder) ;
    for i = 1:2:size(new_mat,2)
        violins([i]).ViolinColor(:) = [1 0 0];
        violins([i+1]).ViolinColor(:) = [0 0 1];
    end

    ylim([-0.5 max(max(new_mat(:,:)))+2])
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    ax.XLabel.Position = [6.5, -28.4, -1];
    
    title(sprintf('ERs evoked per stimulation pair'),'FontSize', 15, 'FontWeight', 'bold')
    ylabel('Number of ERs evoked per stimulation pair','FontSize', 15, 'FontWeight', 'bold')

    legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 12, 'FontWeight', 'bold')

    % Save figure
    outlabel=sprintf('ERsPerStimp_violin.jpg');
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')
    
end
    