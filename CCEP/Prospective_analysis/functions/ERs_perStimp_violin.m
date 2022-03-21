function ERs_perStimp_violin(dataBase,myDataPath)
% Create a violin plot for the number of ERs evoked per stimulation pair. 

% Pre-allocation of matrix with the number of ERs of every session of every
% patient.
new_mat = NaN([54,12]);
subjects = cell(1,size(dataBase,2));
    
for subj = 1:size(dataBase,2)
    Clin = dataBase(subj).agreement_parameter.ERs_stimpClin;
    Prop = dataBase(subj).agreement_parameter.ERs_stimpProp;   
        
     i = 1;
     clin_colm = 2*subj-1;                                  % prealloction of the column number
     prop_colm = 2*subj;                                    % prealloction of the column number

        for stimp = 1:size(Clin,1)                          % For each stimpair                                          
            new_mat(i,clin_colm) = Clin(stimp);             % write the SPES-clin in column 1
            new_mat(i,prop_colm) = Prop(stimp);             % write the SPES-prop in column 2
            i = i+1;           
        end  

        subjects{subj} = dataBase(subj).ccep_clin.sub_label;
end
    
% Create a violin plot
    figure('Position',[205,424,1530,638]);
    violins = violinplot(new_mat) ;                    
    for i = 1:2:size(new_mat,2)
        violins(i).ViolinColor(:) = [1 0 0];
        violins(i+1).ViolinColor(:) = [0 0 1];
    end

    
 % Set significance in plot with a line and not in name
 count = 1;
 ymax = max(max(new_mat));
 
for subj=1:size(dataBase,2)
    if dataBase(subj).statistics.p_ERsperStimp   < 0.01 
        text(count+0.5,ymax-0.2,'**','FontSize',20,'FontWeight','bold')
        plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)


    elseif dataBase(subj).statistics.p_ERsperStimp  < 0.05 
        text(count+0.5,ymax-0.2,'*','FontSize',20,'FontWeight','bold')
        plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)

    end
    count = count+2;
end

     
    % Set labels and axis to own preferences
    ylim([-0.5 max(max(new_mat(:,:)))+2])
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    ax.XTick = 1.5:2:size(new_mat,2);
    ax.XTickLabel = subjects;

    title(sprintf('ERs evoked per stimulation pair'),'FontSize', 15, 'FontWeight', 'bold')
    ylabel('ERs evoked per stimulation pair','FontSize', 15, 'FontWeight', 'bold')

    legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 12, 'FontWeight', 'bold','Position',[0.78,0.80,0.12,0.07])

    medians = median(new_mat,'omitnan');
    ymin = min(ylim);
    y_range = diff(ylim);
    x_as = 1:size(new_mat,2);
    size_pat = size(new_mat,2); 
    second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
    text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')


    % Save figure
    outlabel=sprintf('ERsPerStimp_violin.jpg');
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')
    
end
    