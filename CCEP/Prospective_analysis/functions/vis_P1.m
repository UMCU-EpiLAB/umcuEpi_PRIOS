function vis_P1(myDataPath,dataBase)
% script to determine the statistical value and to visualise the P1 latency

%% Statistics

% Load the presaved P1_latency. Saved during
% PROS01_pipeline_preproces.m
for subj = 1:size(dataBase,2)
    sub_name = dataBase(subj).sub_label;
    P1(subj) = load(fullfile([myDataPath.CCEPpath,'Visualise_agreement/N1_compare/P1_latency/'],[sub_name,'_P1_latency.mat']));
end

% Convert sample to millisecond
fs = 2048;
ts = 1/fs;

for subj = 1:size(dataBase,2)
   
    P1(subj).latency_P1.SPES_clin = P1(subj).latency_P1.SPES_clin* ts *1000 ;
    P1(subj).latency_P1.SPES_prop = P1(subj).latency_P1.SPES_prop* ts *1000 ;
    P1(subj).latency_P1.percentiles_P1_clin = P1(subj).latency_P1.percentiles_P1_clin* ts *1000 ;
    P1(subj).latency_P1.percentiles_P1_prop = P1(subj).latency_P1.percentiles_P1_prop* ts *1000 ;
      
end


for subj = 1:size(dataBase,2)
    %preallocate
    lat_P1 = [];
    
    subName = extractAfter(dataBase(subj).sub_label,'sub-');
    % Make a new 2 column double with the clin value in col 1 and prop
    % value in col 2
    
    lat_P1(:,1) = P1(subj).latency_P1.SPES_clin;
    lat_P1(:,2) = P1(subj).latency_P1.SPES_prop;
    
    % Paired, non-prametric, comparison of two measurements
    % Determine whether the median of the two measurements are
    % significantly different. 
    p_P1 = signrank(lat_P1(:,1), lat_P1(:,2));
    
    dataBase(subj).ccep_clin.p_P1 = p_P1;    
    dataBase(subj).ccep_prop.p_P1 = p_P1;

    % Display the p value 
    fprintf('P1 gives p = %1.4f, for %s. \n', p_P1,subName);
end


%% Visualise in violin plots
% Preallocation
lat_p1 = NaN(457,12); 

for subj = 1:size(dataBase,2)
    
    % Preallocate clinical and propofol columns
    clin_colm = 2*subj-1;                      % prealloction of the column number
    prop_colm = 2*subj;                        % prealloction of the column number

    % create matrix with fall and rise times of clinical and propofol 
    lat_p1(1:size(P1(subj).latency_P1.SPES_clin,1), clin_colm) = P1(subj).latency_P1.SPES_clin;
    lat_p1(1:size(P1(subj).latency_P1.SPES_prop,1), prop_colm) = P1(subj).latency_P1.SPES_prop;
end

% Create violin plot
% make array with medians of fall and rise time with clin als prop
% alternating
medians_P1 = zeros(1,12);  

% create array with medians of latency P1 and latency P0
for subj = 1:size(dataBase,2)
    median_p1_clin(1,subj) = P1(subj).latency_P1.percentiles_P1_clin(2);  
    median_p1_prop(1,subj) = P1(subj).latency_P1.percentiles_P1_prop(2);      
end

medians_P1(1:2:end) = median_p1_clin;
medians_P1(2:2:end) = median_p1_prop;

median_clin = prctile(median_p1_clin,[25 50 75]);
median_prop = prctile(median_p1_prop(2:5), [25 50 75]);
fprintf('median clin P1 = %1.1f ms, percentile = [%1.1f - %1.1f]\n',median_clin(2),median_clin(1),median_clin(3));
fprintf('median prop P1 = %1.1f ms, percentile = [%1.1f - %1.1f]\n',median_prop(2),median_prop(1),median_prop(3));


% Make Violin plot
mat = lat_p1;
medians = medians_P1;

figure('Position',[205,424,1530,638]);

violins = violinplot(mat) ;
for i = 1:2:size(mat,2)
    violins(1,i).ViolinColor = [1 0 0];
    violins(1,i+1).ViolinColor = [0 0 1];
end
    
% Set significance in plot with a line and not in name
count = 1;
ymax = max(max(mat));

 for subj=1:size(dataBase,2)
      p = dataBase(subj).ccep_clin.p_P1;

    if p   < 0.01 
        text(count+0.5,ymax-0.2,'**','FontSize',20,'FontWeight','bold')
        plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)

    elseif p  < 0.05 
        text(count+0.5,ymax-0.2,'*','FontSize',20,'FontWeight','bold')
        plot(count+0.1:0.1:count+0.9, ymax-0.43*ones(9,1),'k','LineWidth',2)

    end
    count = count+2;
 end

 
% Set figure settings
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

% Set double xlabel
ax.XTick = 1.5:2:size(mat,2);
ax.XTickLabel = extractAfter({dataBase(:).sub_label},'sub-');

ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(mat,2);
size_pat = size(mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')

% Draw lines between patients 
ymax = max(ylim); 
for i = 1:2:size(x_as,2)
    x1 = x_as(i)-0.5; 
    if x1 > 0.5
        hold on
        line([x1,x1],[ymin,ymax],'color',[0.8 0.8 0.8]);
    end
end
    
title(sprintf('P1-Latency'),'FontSize', 15, 'FontWeight', 'bold')
ylabel('Latency (ms)','FontSize', 15, 'FontWeight', 'bold')
        
       
legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 12, 'FontWeight', 'bold')
  
% Save figure
outlabel=sprintf('P1_latency_violin.jpg');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/P1_latency/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg') 





end