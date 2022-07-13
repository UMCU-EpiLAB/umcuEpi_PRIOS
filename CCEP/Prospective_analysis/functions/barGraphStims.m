function barGraphStims(dataBase,myDataPath)

%% Make horizontal bar graph for the number of CCEPs per stimulation piar
% Take the ranking of the stimulation pairs based on the clinincal-SPES
% separate per patient
start_row = 1;

for subj = 1:size(dataBase,2)

    if dataBase(subj).ccep_clin.Ckappa <0.6 || dataBase(subj).ccep_prop.Ckappa < 0.6
        % Skip because inter observer agreement is too low

    else
        lastrow = size(dataBase(subj).agreement_parameter.ERs_stimpClin,1) + start_row-1;

        % Clin is left, prop is right
        dataBarPlot(start_row:lastrow, 2) = dataBase(subj).agreement_parameter.ERs_stimpClin *(-1);   
        dataBarPlot(start_row:lastrow, 1) = dataBase(subj).agreement_parameter.ERs_stimpProp;
        
        % use patlabel in the third column to later use to plot in the bar
        % plots when all patients are combined
        dataBarPlot(start_row:lastrow, 3) = str2double(regexp(dataBase(subj).ccep_clin.sub_label,'\d*','Match'));
        
        start_row = size(dataBarPlot,1)+1;
       
    end
end
        
        % Sort based on the clinical-SPES results in the second column. 
        dataBarPlot_sorted = sortrows(dataBarPlot,2,'descend');

        figure('Position',[719,2,1201,1055])
        b1 = barh(dataBarPlot_sorted(:,2),'FaceColor',[17/255, 145/255, 250/255],'EdgeColor',[17/255, 145/255, 250/255]);          % Left --> Clinical
        
        hold on
        b2 = barh(dataBarPlot_sorted(:,1), 'FaceColor',[252/255, 96/255, 57/255],'EdgeColor',[252/255, 96/255, 57/255]);            % Right --> propofol
        
        % legend('prop','clin','Location','southeast')
        
        
%         Now alter the ticks.
        xticks = get(gca, 'xtick');
        
        % Get the current labels
        labels = get(gca, 'xtickLabel'); 
        
        if ischar(labels)
            labels = cellstr(labels); 
        end
        
        % Figure out which ones we need to change
        toscale = xticks < 0;
        
        % Replace the text for the ones < 0
        labels(toscale) = arrayfun(@(x)sprintf('%1.0f', x), ...
                                   abs(xticks(toscale) * -1 ), 'uniformoutput', false);
        
        % Update the tick locations and the labels
        set(gca, 'xtick', xticks, 'xticklabel', labels)
        
        xmax_r = max(get(gca, 'xlim'));
        xmax_l = min(get(gca,'xlim'));
        label(1) = text(xmax_r/2, -11.5, 'Propofol SPES','FontSize',13,'FontWeight','bold');
        label(2) = text(xmax_l/2, -11.5, 'Clinical SPES','FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
        
        title(sprintf('Number of CCEPs detected per stimulation pair for all patients'))
        xlabel('Number of CCEPs detected','HorizontalAlignment','left')
        set(gca,'ycolor','none','XGrid','on','Box','off')
        
        %% Make little txt note in each bar that indicates the number of CCEPs
%         xtips1 = b1.YEndPoints -1.3;
%         ytips1 = b1.XEndPoints;
%         labels1 = arrayfun(@(x)sprintf('%1.0f', x), ...
%                                    abs(b1.YData * -1 ), 'uniformoutput', false);
%         text(xtips1,ytips1,labels1,'VerticalAlignment','middle')
%         
%         xtips2 = b2.YEndPoints +0.3;
%         ytips2 = b2.XEndPoints;
%         labels2 = string(b2.YData);
%         text(xtips2,ytips2,labels2,'VerticalAlignment','middle')
% 
% 
%         % Subject label in the bars
%         xtipssubj = ones(1,size(xtips1,2))*-1;              % just a little left of the midline.
%         subLabels1 = strings(size(xtips1,2),1);
%         for r =1:size(xtips1,2)
%             subLabels1(r,1) = sprintf('0%1.0f', dataBarPlot_sorted(r,3));          
%         end
% 
%         text(xtipssubj,ytips1,subLabels1,'VerticalAlignment','middle','Color',[1 1 1],'FontSize',4)
%         text(xtipssubj*-1,ytips1,subLabels1,'VerticalAlignment','middle','Color',[1 1 1],'FontSize',4)
% 
%     


% Save figure
outlabel=sprintf('allsubs_barGraph.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/BarGraph/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')


%% Bar graph for indegree and outdegree

%% Make horizontal bar graph for the number of CCEPs per stimulation piar
% Take the ranking of the stimulation pairs based on the clinincal-SPES
% separate per patient

mode = {'In-degree','Out-degree','Betweenness Centrality'};
figure('Position',[719,2,1201,1055])


for m = 1:size(mode,2)
    
    start_row = 1;
    dataBarPlot = [];
    dataBarPlot_sorted = [];
    
    for subj = 1:size(dataBase,2)

        if isequal(mode{m},'In-degree')
            data_clin = dataBase(subj).agreement_parameter.indegreeN_Clin;
            data_prop = dataBase(subj).agreement_parameter.indegreeN_Prop;

        elseif isequal(mode{m},'Out-degree')
            data_clin = dataBase(subj).agreement_parameter.outdegreeN_Clin;
            data_prop = dataBase(subj).agreement_parameter.outdegreeN_Prop;
    
        elseif isequal(mode{m},'Betweenness Centrality')
            data_clin = dataBase(subj).agreement_parameter.BCN_Clin;
            data_prop = dataBase(subj).agreement_parameter.BCN_Prop;

        end


    
        if dataBase(subj).ccep_clin.Ckappa <0.6 || dataBase(subj).ccep_prop.Ckappa < 0.6
            % Skip because inter observer agreement is too low
    
        else
            lastrow = size(data_clin,2) + start_row-1;
    
            % Clin is left, prop is right
            dataBarPlot(start_row:lastrow, 2) = data_clin *(-1);   
            dataBarPlot(start_row:lastrow, 1) = data_prop;
            
            % use patlabel in the third column to later use to plot in the bar
            % plots when all patients are combined
            dataBarPlot(start_row:lastrow, 3) = str2double(regexp(dataBase(subj).ccep_clin.sub_label,'\d*','Match'));
            
            start_row = size(dataBarPlot,1)+1;
           
        end
    end
            
            % Sort based on the clinical-SPES results in the second column. 
            dataBarPlot_sorted = sortrows(dataBarPlot,2,'descend');
    
            subplot(3,1,m)
            b1 = barh(dataBarPlot_sorted(:,2),'FaceColor',[17/255,145/255,250/255],'EdgeColor',[17/255,145/255,250/255]);     % Left --> Clinical
            
            hold on
            b2 = barh(dataBarPlot_sorted(:,1), 'FaceColor',[252/255,96/255,57/255],'EdgeColor',[252/255,96/255,57/255]);    % Right --> propofol
            
            % legend('prop','clin','Location','southeast')
            
            
    %         Now alter the ticks.
            xticks = get(gca, 'xtick');
          
            % Get the current labels
            labels = get(gca, 'xtickLabel'); 
            
            if ischar(labels)
                labels = cellstr(labels); 
            end
            
            % Figure out which ones we need to change
            toscale = xticks < 0;
            
            % Replace the text for the ones < 0
            labels(toscale) = arrayfun(@(x)sprintf('%1.2f', x), ...
                                       abs(xticks(toscale) * -1 ), 'uniformoutput', false);
            
            % Update the tick locations and the labels
            set(gca, 'xtick', xticks, 'xticklabel', labels)

%             if m == 3
% 
%                 xmax_r = max(get(gca, 'xlim'));
%                 xmax_l = min(get(gca,'xlim'));
%                 label(1) = text(xmax_r/2, -55.5, 'Propofol SPES','FontSize',13,'FontWeight','bold');
%                 label(2) = text(xmax_l/2, -55.5, 'Clinical SPES','FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
%             end
               
            legend({'Clinical-SPES','Propofol-SPES'},'FontSize',8,'Location','southeast','Box','off')

            % To determine the significance for the ranking
            % Use Spearman 'rows','pairwise' to ensure that row's with NaN's in both columns are not considered in the analysis.
            [rho, pval] = corr(dataBarPlot(:,1), -dataBarPlot(:,2),'Type','Spearman','rows','pairwise');      
            
            if pval < 0.001
                    title(sprintf('%s, p = <0.001, r_s = %1.3f', mode{m}, rho),'FontSize',12)               
            else
                    title(sprintf('%s, p = %1.3f, r_s = %1.3f', mode{m}, pval, rho),'FontSize',12)
            end 

%             xlabel(sprintf('Normalised %s',mode{m}),'HorizontalAlignment','left')
            set(gca,'ycolor','none','XGrid','on','Box','off')

            
end

% Save figure
outlabel=sprintf('allsubs_barGraph.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/BarGraph/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

end