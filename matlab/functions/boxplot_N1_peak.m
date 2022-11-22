function [dataBase] = boxplot_N1_peak(dataBase, myDataPath)
% Make boxplots of the latency of the N1 peaks
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

%% Make VIOLINPlots    
% New_mat is a table that contains all latencies per patient. This is NOT
% electrode specific. It is just a concatenation of all N1's detected
new_mat = NaN(451,size(dataBase,2)*2);  
pval = NaN(size(dataBase,2,1));

fs = 2048;
ts = 1/fs; 

for subj = 1:size(dataBase,2)

    clin = dataBase(subj).ccep_clin.n1_peak_sample;
    clin = ((clin*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
    prop = dataBase(subj).ccep_prop.n1_peak_sample;
    prop = ((prop*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact

    dataBase(subj).lat_elec_clin = [];
    dataBase(subj).lat_elec_prop = [];

    i = 1;     

    for stimp = 1:size(dataBase(subj).ccep_prop.stimsets_avg,1)                          % For each stimpair
        for elec = 1:size(dataBase(subj).ccep_prop.ch,1)                                   % For each electrode

        % When both clinical SPES AND propofol SPES show an ER
          if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
              % Save latency of each response that occurs in clin and prop
              % per patient. This is used to determine the correlation 
              % between the latency. DEPTH ELECTRODES in response and stimpair are excluded.
              elec1_stimp = dataBase(subj).ccep_prop.stimsets_avg(stimp,1);
              elec2_stimp = dataBase(subj).ccep_prop.stimsets_avg(stimp,2);

              if ~isequal(dataBase(subj).ccep_clin.tb_channels.group{elec,:}, 'depth') && ~isequal(dataBase(subj).ccep_clin.tb_channels.group{elec1_stimp,:}, 'depth') && ~isequal(dataBase(subj).ccep_clin.tb_channels.group{elec2_stimp,:}, 'depth') 
                  dataBase(subj).lat_elec_clin(i,1) = clin(elec, stimp);
                  dataBase(subj).lat_elec_prop(i,1) = prop(elec, stimp);
                  i = i+1;
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

% Replace zero with NaN to avoid influence on the mean
new_mat((new_mat == 0)) = NaN;                                      

% Add 2 extra columns with the collection of all responses of ALL PATIENTS COMBINED during
% clinical-SPES and during propofol-SPES
new_mat(1:2706,13) = [new_mat(:,1);new_mat(:,3);new_mat(:,5);new_mat(:,7);new_mat(:,9);new_mat(:,11)];
new_mat(1:16236,14) = [new_mat(:,2);new_mat(:,4);new_mat(:,6);new_mat(:,8);new_mat(:,10);new_mat(:,12)];

new_mat(new_mat == 0)= NaN;

% Create boxplot with the amplitude of SPES clin and SPES prop
figure('Position',[205,424,1530,638]);

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

% Set double xlabel
ax.XTick = 1.5:2:size(new_mat,2);
pat_names = arrayfun(@(x) x.ccep_clin.sub_label, dataBase, 'UniformOutput', false);
ax.XTickLabel  = pat_names'; 
ax.XTickLabel(end+1) = {'Patients combined'};

medians = median(new_mat,'omitnan');

ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(new_mat,2);
size_pat = size(new_mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.1f '),' '));
text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12)

% Draw lines between separate patient results and all patients combined
ymax = max(ylim); 
x1 = x_as(13)-0.5; 
hold on
line([x1,x1],[ymin,ymax],'color','k','LineWidth',3,'LineStyle','--');

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

end
    