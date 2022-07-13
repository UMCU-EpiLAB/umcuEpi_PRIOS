function [dataBase] = boxplot_N1_peak(dataBase, myDataPath)
% Make boxplots of the latency of the N1 peaks.


% Remove  subjects with too low interobserver agreement
dataBase_remove = zeros(1,size(dataBase,2));
for s = 1:size(dataBase,2)
     if dataBase(s).ccep_clin.Ckappa <0.6 || dataBase(s).ccep_prop.Ckappa < 0.6
            % Skip because inter observer agreement is too low
        dataBase_remove(:,s) = 1;
     else
        dataBase_remove(:,s) = 0;
     end
end

dataBase(dataBase_remove ==1) = [];
fs = 2048;
ts = 1/fs; 


% % Statistics
% for subj = 1:size(dataBase,2)
%    
%     ccep_clin = dataBase(subj).ccep_clin;
%     ccep_prop = dataBase(subj).ccep_prop;
% 
%     %ts = 1/(size(dataBase_prop.tt,2)/4);                             % Devide tt by four because ccep_prop.tt includes 4 seconds.
%     fs = 2048;
%     ts = 1/fs; 
% 
%     clin = ccep_clin.n1_peak_sample;
%     clin = ((clin*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
%     prop = ccep_prop.n1_peak_sample;
%     prop = ((prop*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
% 
% 
%     % Create matrix with the clinical values of the latency in the first column and
%     % the propofol values in the second
%     i = 1;
%     for stimp = 1:size(ccep_prop.stimsets_avg,1)                          % For each stimpair
%         for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode
% 
%         % When both clinical SPES and propofol SPES show an ER
%           if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
%                 new_mat(i,1) = clin(elec, stimp);            %#ok<AGROW> % plot the SPES-clin amp in column 1
%                 new_mat(i,2) = prop(elec, stimp);            %#ok<AGROW> % plot the SPES-prop amp in column 2
%                 i = i+1;                
%           end
%         end      
%     end
% 
%     % Determine gaussianity
%     % Biggest part of the data is not normally distributed (determined on 10-12-2020)
% %         NorDisClin = lillietest(new_mat(:,1))   ;               % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected 
% %         NorDisProp = lillietest(new_mat(:,2));
% 
%     % Paired, non-prametric, comparison of two measurements
%     % Determine whether the median of the two measurements are
%     % significantly different. 
%     p = signrank(new_mat(:,1), new_mat(:,2)) ;                       
%     dataBase(subj).ccep_clin.p_n1 = p;
%     dataBase(subj).ccep_clin.mean_N1_lat = mean(new_mat(:,1));
%     dataBase(subj).ccep_prop.p_n1 = p;
%     dataBase(subj).ccep_prop.mean_N1_lat = mean(new_mat(:,2));
% 
%     % Display the p value 
%     if p<0.05
%         fprintf('Test between the N1-Latency of %s SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is a significant difference between the two protocols \n',ccep_clin.sub_label, p);
%         fprintf('Median clin = %1.2f, Median prop = %1.2f\n', dataBase(subj).ccep_clin.mean_N1_lat ,dataBase(subj).ccep_prop.mean_N1_lat)
%     else
%         fprintf('Test between the N1-Latency of %s SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is NO significant difference between the two protocols \n',ccep_clin.sub_label, p);
%         fprintf('Median clin = %1.2f, Median prop = %1.2f \n', dataBase(subj).ccep_clin.mean_N1_lat ,dataBase(subj).ccep_prop.mean_N1_lat)
% 
%     end
% end


%% Make VIOLINPlots    
% New_mat is a table that contains all latencies per patient. This is NOT
% electrode specific. It is just a concatenation of all N1's detected
new_mat = [];  

for subj = 1:size(dataBase,2)

    ccep_clin = dataBase(subj).ccep_clin;
    ccep_prop = dataBase(subj).ccep_prop;

    clin = ccep_clin.n1_peak_sample;
    clin = ((clin*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
    prop = ccep_prop.n1_peak_sample;
    prop = ((prop*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact

    dataBase(subj).lat_elec_clin = [];
    dataBase(subj).lat_elec_prop = [];

     i = 1;     

    for stimp = 1:size(ccep_prop.stimsets_avg,1)                          % For each stimpair
        for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode

        % When both clinical SPES AND propofol SPES show an ER
          if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
%               new_mat(i,clin_colm) = clin(elec, stimp);            % plot the SPES-clin amp in column 1
%               new_mat(i,prop_colm) = prop(elec, stimp);            % plot the SPES-prop amp in column 2          

              % Save latency of each response that occurs in clin and prop
              % per patient. This is used to determine the correlation 
              % between the distance and latency. Therefore DEPTH
              % ELECTRODES in response and stimpair are excluded.
              elec1_stimp = ccep_prop.stimsets_avg(stimp,1);
              elec2_stimp = ccep_prop.stimsets_avg(stimp,2);

              if ~isequal(ccep_clin.tb_channels.group{elec,:}, 'depth') && ~isequal(ccep_clin.tb_channels.group{elec1_stimp,:}, 'depth') && ~isequal(ccep_clin.tb_channels.group{elec2_stimp,:}, 'depth') 
                  dataBase(subj).lat_elec_clin(i,1) = clin(elec, stimp);
                  dataBase(subj).lat_elec_prop(i,1) = prop(elec, stimp);
                  i = i+1;
              else
%                   dataBase(subj).lat_elec_clin(i,1) = NaN;
%                   dataBase(subj).lat_elec_prop(i,1) = NaN;
%                   i = i+1;
              end

               
          end        

        end      
    end 

    % Wilcoxon signed rank test
    pval(subj,:) = signrank(dataBase(subj).lat_elec_clin(:,1) , dataBase(subj).lat_elec_prop(:,1));

    % create new_mat with latencies of all patients used for plotting
    clin_colm = 2*subj-1;                      % prealloction of the column number
    prop_colm = 2*subj;                        % prealloction of the column number
    new_mat(1:size(dataBase(subj).lat_elec_clin(:,1),1), clin_colm) = dataBase(subj).lat_elec_clin(:,1);
    new_mat(1:size(dataBase(subj).lat_elec_clin(:,1),1), prop_colm) = dataBase(subj).lat_elec_prop(:,1);

end



new_mat((new_mat == 0)) = NaN;                                      % replace zero with NaN to avoid influence on the mean
% av_lat_elec((av_lat_elec== 0)) = NaN;                                      % replace zero with NaN to avoid influence on the mean


% %% Make table of new_mat
% % For easier interpretation of the matrix and this table can be saved and
% % later used in figures where latency is plotted on the brain
% 
% % Perallocation
% sz = [size(new_mat,1) size(dataBase,2)*2];
% varTypes = repmat({'double'},size(dataBase,2)*2,1);
% table_latency = table('Size',sz,'VariableTypes',varTypes);
% pat_names = cell(size(dataBase,2),1);
% 
% % Make an array with all patient names to use in the figures and table
% for pat = 1:size(dataBase,2)
%     pat_names{pat,:} = dataBase(pat).ccep_clin.sub_label; 
% end
% 
% % Give each column the right name (patient and clin/prop)
% for col = 1:size(pat_names,1)
%     clin_colm = 2*col-1;                      
%     prop_colm = 2*col; 
%     % Add column name
%     table_latency.Properties.VariableNames{clin_colm} = [pat_names{col},'_clin'];
%     table_latency.Properties.VariableNames{prop_colm} = [pat_names{col},'_prop'];
%     % Add latency to table 
%     table_latency(:,clin_colm) = table(new_mat(:,clin_colm));
%     table_latency(:,prop_colm) = table(new_mat(:,prop_colm));
% end
% 

%%
% medians =  median(new_mat,'omitnan');

% Create boxplot with the amplitude of SPES clin and SPES prop
figure('Position',[205,424,1530,638]);
% columnMeans = mean(new_mat, 1, 'omitnan');

violins = violinplot(new_mat) ;
for i = 1:2:size(new_mat,2)
    violins(1,i).ViolinColor = [17/255, 145/255, 250/255];
    violins(1,i+1).ViolinColor = [252/255, 96/255, 57/255];
end

% Set significance in plot with a line and not in name
count = 1;
ymax = max(max(new_mat));
 
for subj=1:size(dataBase,2)

        if pval(subj)   < 0.001 
            text(count+0.5,ymax-0.2,'**','FontSize',20)
            plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)
        
        
        elseif pval(subj) < 0.05 
            text(count+0.5,ymax-0.2,'*','FontSize',20)
            plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)
                
        end
        count = count+2;
end


ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
%ax.XAxis.FontWeight = 'bold';
%ax.YAxis.FontWeight = 'bold';

% Set double xlabel
ax.XTick = 1.5:2:size(new_mat,2);
pat_names = arrayfun(@(x) x.ccep_clin.sub_label, dataBase, 'UniformOutput', false);
ax.XTickLabel  = pat_names'; 

medians = median(new_mat,'omitnan');

ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(new_mat,2);
size_pat = size(new_mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12)

% Draw lines between patients 
ymax = max(ylim); 
for i = 1:2:size(x_as,2)
    x1 = x_as(i)-0.5; 
    if x1 > 0.5
        hold on
        line([x1,x1],[ymin,ymax],'color',[0.8 0.8 0.8]);
    end
end

title(sprintf('N1 Latency'),'FontSize', 15)
ylabel('Latency (milliseconds)','FontSize', 15)
legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 10,'Position',[0.804052288487456 0.796238244514106 0.0941176452473098 0.072194357366771])

   
% Save figure
outlabel=sprintf('Latency_violin.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')


%% grouping all patients together and determine the Wilcoxon singed rank test
% To compare the medians of all patients combined
% Concatenate the N1-peak latencies of all responses during clinical spes.
    % Remember that only responses are considered that showed an N1 in
    % clinical and propofol spes, therefore the arrays of clin and prop
    % should be of equal size.
all_pat_clin = new_mat(:,[1,3,5,7,9,11]);
all_pat_clin = all_pat_clin(:);
all_pat_clin(isnan(all_pat_clin)) = [];

all_pat_prop = new_mat(:,[2,4,6,8,10,12]);
all_pat_prop = all_pat_prop(:);
all_pat_prop(isnan(all_pat_prop)) = [];

% Determine comparison off N1-peak latencies for all patients combined
pval_allPat = signrank(all_pat_clin , all_pat_prop);

prctile(all_pat_clin,[25 50 75])
prctile(all_pat_prop,[25 50 75])



%%
% % % % %% Calculate the median N1-peak latency per electrode and then compare the N1-latency change per electrode during clin and prop
% % % % 
% % % % for subj = 1:size(dataBase,2)
% % % % 
% % % %     med_per_elec = [];
% % % %     ccep_clin = dataBase(subj).ccep_clin;
% % % %     ccep_prop = dataBase(subj).ccep_prop;
% % % % 
% % % %     clin = ccep_clin.n1_peak_sample;
% % % %     clin = ((clin*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
% % % %     prop = ccep_prop.n1_peak_sample;
% % % %     prop = ((prop*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
% % % % 
% % % % %     clin_colm = 2*subj-1;                      % prealloction of the column number
% % % % %     prop_colm = 2*subj;                        % prealloction of the column number
% % % % 
% % % %     % determine the median for all reponses per electrode 
% % % %     % independent of whether the same stimulation pairs evoked the reponses
% % % %     for elec = 1:size(clin,1)                                   % For each electrode
% % % %           
% % % %         med_per_elec(elec,1) = median(clin(elec,:),'omitnan');            % plot the SPES-clin amp in column 1
% % % %         med_per_elec(elec,2) = median(prop(elec,:),'omitnan');            % plot the SPES-prop amp in column 2               
% % % % 
% % % %     end
% % % % 
% % % % 
% % % %     [rho, pval] = corr(med_per_elec(:,1), med_per_elec(:,2),'Type','Spearman','rows','pairwise');      
% % % % 
% % % %     fprintf('PRIOS0%1.0f, p = %1.3f \n', subj, pval)
% % % % 
% % % % end



        
%% Make scatter of the latency
% fig = figure('Position',[302,17,938,1039]);
% 
% for subj = 1:size(dataBase,2)
%         
%        colm_clin = 2*subj-1;
%        colm_prop = (2*subj);
%        clin = new_mat(:,colm_clin);
%        prop = new_mat(:,colm_prop);
%        p = dataBase(subj).ccep_clin.p_n1;
%        
%        % do not plot values with NaN
%        if size(clin(isnan(clin)),1) > 1
%           clin = clin(~isnan(clin));
%           prop = prop(~isnan(prop));              
%        end
%        
%         subplot(size(dataBase,2),1,subj)
%         plot(clin  , prop,'*'); 
%         title(sprintf('%s, p =  %0.3f', dataBase(subj).ccep_clin.sub_label, p),'FontSize',12);
%         % Change fontsize
%         ax = gca;
%         ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 14;
%         if p < 0.01
%              title(sprintf('%s, p = <0.01', dataBase(subj).ccep_clin.sub_label),'FontSize',12)
%          elseif p<0.05
%              title(sprintf('%s, p = <0.05', dataBase(subj).ccep_clin.sub_label),'FontSize',12)
%         end 
%                  
%   
% %             if p< 0.05
% %                 
% %                 [P,S] = polyfit(clin,prop,1);
% %                 [y_fit, ~] = polyval(P,clin,S);
% % 
% %                 plot(clin,prop,'*')                        % This is equal to scatter
% %                 hold on
% % 
% %                 % Plot polyfit throught data points
% %                 plot(clin,y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
% %                 hold on
% %                 
% %                  
% %             
% %             end
%            
%         % Xmax is determined by the value on the x-axis which is the
%         % SPES-clin
%         xmin = min(clin);                 
%         xmax = max(clin);   
%         xlim([xmin, xmax]) 
%                            
% end 
% 
% % Create main title without decreasing the subplots sizes
% annotation(gcf,'textbox',[0.32 0.95 0.35 0.043],'VerticalAlignment','middle','String',sprintf('N1-Latency'),'HorizontalAlignment','center','FontWeight','bold','FontSize',15,'FitBoxToText','off','EdgeColor','none');   
%      
% %Create one main x- and y-label
% han=axes(fig,'visible','off'); 
% han.YLabel.Visible='on';
% han.YLabel.Position = [-0.0584    0.5000         0];
% ylabel(han,{'Value SPES-prop'}, 'FontSize',17,'fontweight','bold')
% 
% han.XLabel.Visible='on';
% xlabel(han,{'Value SPES-clin'}, 'FontSize',17,'fontweight','bold')
% 
%      
% 
% % Save figure
% outlabel=sprintf('n1_scatter_Latency.png');
% path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/Scatter/');
% if ~exist(path, 'dir')
%     mkdir(path);
% end
% saveas(gcf,[path,outlabel],'png')
        
%% Determine the multiplication factor

% for subj = 1:size(dataBase,2)
%   M_Clin = medians(1:2:end);
%   M_Prop = medians(2:2:end);
% end
% 
% %%%% Remove this (just for the report)
% per_N1_clin = prctile(M_Clin, [25 50 75])               
% per_N1_prop= prctile(M_Prop, [25 50 75])
% 
% 
% % First column contains median N1-latency per patient during clinical-SPES
% % Second column contains median N1-latency per patient during propogol-SPES
% T_N1(:,1) = medians(1:2:end);
% T_N1(:,2) = medians(2:2:end);
% 
% % Third column contains the multiplication factor between propofol-SPES and
% % clinical-SPES (The number of times propofol was larger than clinical)
% for i = 1:size(T_N1,1)
%     T_N1(i,3) = T_N1(i,2)/T_N1(i,1);
% end
% 
% % Last ROW contains the median N1-latency during clinical SPES and propofl
% % SPES for all patients
% Size_mat = (size(T_N1,1)+1);
% T_N1(Size_mat,1) = median(M_Clin);
% T_N1(Size_mat,2) = median(M_Prop);
% T_N1(Size_mat,3) = NaN;                 % not required to calculate a median of the multiplication factor
% % T_N1(Size_mat,3) = median(T_N1(1:size(M_Clin,2),3));
% 
% variables = {'N1 clinical','N1 propofol','Mult-factor'};
% 
% T_N1 = table(T_N1(:,1),T_N1(:,2),T_N1(:,3), 'VariableNames',variables,'RowNames',[pat_names','Medians']);
% disp(T_N1)
end
    