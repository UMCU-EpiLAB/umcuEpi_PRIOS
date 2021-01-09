function ERs_perStimp_violin(dataBase,myDataPath)
% Create a violin plot for the number of ERs evoked per stimulation pair. 

% Pre-allocation of matrix with the number of ERs of every session of every
% patient.
new_mat = NaN([54,12]);
    
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
end
    
% Create a violin plot
    figure('Position',[205,424,1530,638]);
    
% Create the names with astrics for significance
subj_sig = cell(1,size(dataBase,2));
for subj = 1:size(dataBase,2)
    
    subj_name = extractAfter(dataBase(subj).sub_label,'sub-');
    Sig_num = dataBase(subj).statistics.p_ERsperStimp;  
    
    if Sig_num < 0.01
        subj_sig{subj} = sprintf('%s**',subj_name);
        
    elseif Sig_num<0.05
        subj_sig{subj} = sprintf('%s*', subj_name);
        
    else
        subj_sig{subj} = subj_name;
    end  
       
end
    grouporder =  [subj_sig(1),' ',subj_sig(2),' ',subj_sig(3),' ',subj_sig(4),' ',subj_sig(5),' ',subj_sig(6),' '];

    % Plot violin plot
    violins = violinplot(new_mat,grouporder) ;                    
    for i = 1:2:size(new_mat,2)
        violins(i).ViolinColor(:) = [1 0 0];
        violins(i+1).ViolinColor(:) = [0 0 1];
    end

    % Set labels and axis to own preferences
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
    