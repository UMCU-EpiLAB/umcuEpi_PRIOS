function CCEPS_per_lobe_elec(myDataPath,dataBase, Destrieux_label_pat, roi_central, roi_frontal, roi_occipital, roi_temporal, roi_parietal)

fs = 2048;
ts = 1/fs;

% For each electrode of the stimulation pair, determine the lobe, and
% then make a matrix with all the responses that were evoked BY that
% lobe
r_temp = 1;
r_front = 1;
r_cent = 1;
r_pari = 1;

mode = [{'Clin'}, {'Prop'}];
lat_temporal_elec = struct;
lat_parietal_elec = struct;
lat_frontal_elec = struct;
lat_central_elec = struct;

for pat = 1:size(dataBase,2) 

    elec_lobe_label = cell(size(Destrieux_label_pat(:,pat),1),1);
   
    % Rewrite destrieux_label_pat to the lobes instead of the destrieux
    % labels
    elec_lobe_label(ismember(Destrieux_label_pat(:,pat),roi_central),1) = {'Central'};
    elec_lobe_label(ismember(Destrieux_label_pat(:,pat),roi_frontal),1) = {'Frontal'};
    elec_lobe_label(ismember(Destrieux_label_pat(:,pat),roi_occipital),1) = {'Occipital'};
    elec_lobe_label(ismember(Destrieux_label_pat(:,pat),roi_temporal),1) = {'Temporal'};
    elec_lobe_label(ismember(Destrieux_label_pat(:,pat),roi_parietal),1) = {'Parietal'};
    elec_lobe_label(isnan(Destrieux_label_pat(:,pat))) = [];
    
    % Extra check to test whether the same number of chnnels is still
    % vallid
    if ~isequal(length(elec_lobe_label), size(dataBase(pat).ccep_clin.ch,1))
        error('elec_lobe_label is not the same length as the registered channels')
    end

    
    sz_col = 1:size(dataBase(pat).ccep_clin.n1_peak_sample,2);  % Length of columns is different for every patient, this ensures that the rows are concatenated

 
    % For each stimulation pair determine the lobe it is located on and
    % save the evoked CCEPs on that electrode in a matrix per lobe
    for elec = 1:size(dataBase(pat).ccep_clin.ch,1)                  % Does not matter whether this is clin or prop, is equal
                      
        % Location of stimulation pair electrodes
        lobe_elec = elec_lobe_label(elec,1);
                
        for m = 1:size(mode,2)
            if isequal(mode{m},'Clin')
                data = dataBase(pat).ccep_clin;
                
            elseif isequal(mode{m},'Prop')
                data = dataBase(pat).ccep_prop;

            end

             % Save the responses of the electrode when located per lobe            
             if isequal(lobe_elec{:},'Temporal') 
                lat_temporal_elec.(sprintf('%s',mode{m}))(r_temp, sz_col) = (data.n1_peak_sample(elec,:) -(2*fs)) *ts *1000;
                r_temp = r_temp+1;

            elseif isequal(lobe_elec{:},'Frontal')
                lat_frontal_elec.(sprintf('%s',mode{m}))(r_front, sz_col) = (data.n1_peak_sample(elec,:) -(2*fs)) *ts *1000;          
                r_front = r_front+1;

            elseif isequal(lobe_elec{:},'Parietal')
                lat_parietal_elec.(sprintf('%s',mode{m}))(r_pari, sz_col) = (data.n1_peak_sample(elec,:) -(2*fs)) *ts *1000;
                r_pari = r_pari+1;    
                   
            elseif isequal(lobe_elec{:},'Central')
                lat_central_elec.(sprintf('%s',mode{m}))(r_cent, sz_col) = (data.n1_peak_sample(elec,:) -(2*fs)) *ts *1000;    
                r_cent = r_cent+1;               
             end
        end       
    end
end



% Replace zeros with NaN to avoid influence on median
lat_temporal_elec.Clin(ismember(lat_temporal_elec.Clin,0)) = NaN;
lat_frontal_elec.Clin(ismember(lat_frontal_elec.Clin,0)) = NaN;
lat_parietal_elec.Clin(ismember(lat_parietal_elec.Clin,0)) = NaN;
lat_central_elec.Clin(ismember(lat_central_elec.Clin,0)) = NaN;

lat_temporal_elec.Prop(ismember(lat_temporal_elec.Prop,0)) = NaN;
lat_frontal_elec.Prop(ismember(lat_frontal_elec.Prop,0)) = NaN;
lat_parietal_elec.Prop(ismember(lat_parietal_elec.Prop,0)) = NaN;
lat_central_elec.Prop(ismember(lat_central_elec.Prop,0)) = NaN;

% Make array of latencies in matrix
lat_temporal_elec.Clin(isnan(lat_temporal_elec.Clin)) = [];
lat_temporal_elec.Prop(isnan(lat_temporal_elec.Prop)) = [];
lat_frontal_elec.Clin(isnan(lat_frontal_elec.Clin)) = [];
lat_frontal_elec.Prop(isnan(lat_frontal_elec.Prop)) = [];
lat_parietal_elec.Clin(isnan(lat_parietal_elec.Clin)) = [];
lat_parietal_elec.Prop(isnan(lat_parietal_elec.Prop)) = [];
lat_central_elec.Clin(isnan(lat_central_elec.Clin)) = [];
lat_central_elec.Prop(isnan(lat_central_elec.Prop)) = [];


%% %% Display all latencies per lobe in a violinplot
violin_mat = NaN(1000,8);
violin_mat(1:size(lat_temporal_elec.Clin,2), 1) = lat_temporal_elec.Clin';    % Temporal Clin
violin_mat(1:size(lat_temporal_elec.Prop,2), 2) = lat_temporal_elec.Prop';    % Temporal Prop

violin_mat(1:size(lat_frontal_elec.Clin,2), 3) = lat_frontal_elec.Clin';    % Frontal Clin
violin_mat(1:size(lat_frontal_elec.Prop,2), 4) = lat_frontal_elec.Prop';    % Frontal Prop

violin_mat(1:size(lat_parietal_elec.Clin,2), 5) = lat_parietal_elec.Clin';    % Parietal Clin
violin_mat(1:size(lat_parietal_elec.Prop,2), 6) = lat_parietal_elec.Prop';    % Parietal Prop

violin_mat(1:size(lat_central_elec.Clin,2), 7) = lat_central_elec.Clin';    % Central Clin
violin_mat(1:size(lat_central_elec.Prop,2), 8) = lat_central_elec.Prop';    % Central Prop


figure('Position',[205,424,1530,638]);
violins = violinplot(violin_mat) ;                    
for i = 1:2:size(violin_mat,2)
    violins(1,i).ViolinColor = [1 0 0];
    violins(1,i+1).ViolinColor = [0 0 1];
end

% Calculate the number of CCEPs evoked for each lobe
for i = 1:size(violin_mat,2)
    nr_CCEP(i,:) = sum(~isnan(violin_mat(:,i)));
end

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

medians = median(violin_mat,'omitnan');
ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(violin_mat,2);
size_pat = size(violin_mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.1f '),' '));
text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.08*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12)

third_row_txt = cellstr(strsplit(num2str(nr_CCEP','%1.0f '),' '));
text([(x_as(1)-x_as(2))*0.5 x_as], ones(1,size_pat+1)*ymin-0.12*y_range, ['Total CCEPs' third_row_txt],'HorizontalAlignment','center','FontSize', 12)

% Draw lines between patients 
ymax = max(ylim); 
for i = 1:2:size(x_as,2)
    x1 = x_as(i)-0.5; 
    if x1 > 0.5
        hold on
        line([x1,x1],[ymin,ymax],'color',[0.8 0.8 0.8]);
    end
end

ax.XTick = 1.5:2:size(violin_mat,2);
ax.XTickLabel  = [{'Temporal'} {'Frontal'} {'Parietal'} {'Central'}]; 


title(sprintf('N1-Latency evoked on electrodes per lobe '),'FontSize', 15, 'FontWeight', 'bold')
ylabel('N1-latency of response electrodes per lobe','FontSize', 15, 'FontWeight', 'bold')


% STATISTICS
p_temp = ranksum(lat_temporal_elec.Clin, lat_temporal_elec.Prop);
p_fron = ranksum(lat_frontal_elec.Clin, lat_frontal_elec.Prop);
p_pari = ranksum(lat_parietal_elec.Clin, lat_parietal_elec.Prop);
p_cent = ranksum(lat_central_elec.Clin, lat_central_elec.Prop);


fprintf('P-value between clin and prop Temporal for all patients = %1.3f \n',p_temp);
fprintf('P-value between clin and prop Frontal for all patients = %1.3f \n',p_fron);
fprintf('P-value between clin and prop Parietal for all patients = %1.3f \n',p_pari);
fprintf('P-value between clin and prop Central for all patients = %1.3f \n',p_cent);


% Display lines with significance in plot
ymax = max(max(violin_mat));
count = 1;
if p_temp  < 0.01 
 text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.01:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)


elseif p_temp  < 0.05 
 text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)

end

count = 3;
if p_fron  < 0.01 
 text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)


elseif p_fron  < 0.05 
 text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)

end

count =5 ;
if p_pari  < 0.01 
 text(count+0.5,ymax-0.1,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-0.93*ones(9,1),'k','LineWidth',2)


elseif p_pari  < 0.05 
 text(count+0.5,ymax-0.1,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-0.93*ones(9,1),'k','LineWidth',2)

end

count =7 ;
if p_cent  < 0.01 
 text(count+0.5,ymax-0.1,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-0.93*ones(9,1),'k','LineWidth',2)


elseif p_cent  < 0.05 
 text(count+0.5,ymax-0.1,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-0.93*ones(9,1),'k','LineWidth',2)

end

legend([violins(1).ViolinPlot,violins(2).ViolinPlot], 'Clinical SPES','Propofol SPES','FontSize', 12, 'FontWeight', 'bold','Position',[0.78,0.845,0.12,0.07])

    
    % Save figure
outlabel=sprintf('N1-latency_lobe_elec.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')





end
