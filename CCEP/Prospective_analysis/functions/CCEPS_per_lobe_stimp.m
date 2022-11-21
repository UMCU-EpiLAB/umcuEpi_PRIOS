function CCEPS_per_lobe_stimp(myDataPath,dataBase, Destrieux_label_pat, roi_central, roi_frontal, roi_occipital, roi_temporal, roi_parietal)

fs = 2048;
ts = 1/fs;

% For each electrode of the stimulation pair, determine the lobe, and
% then make a matrix with all the responses that were evoked BY that
% lobe
c_temp = 1;
c_front = 1;
c_cent = 1;
c_pari = 1;

mode = [{'Clin'}, {'Prop'}];
lat_temporal_stimp = struct;
lat_parietal_stimp = struct;
lat_frontal_stimp = struct;
lat_central_stimp = struct;

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

    
    sz_row = 1:size(dataBase(pat).ccep_clin.n1_peak_sample,1);  % Length of rows is different for every patient, this ensures that the columns are concatenated

 
    % For each stimulation pair determine the lobe it is located on and
    % save the evoked CCEPs by that stimulation pair in a matrix per lobe
    for stimp = 1:size(dataBase(pat).ccep_clin.stimsets_avg,1)                  % Does not matter whether this is clin or prop, is equal
        stimp1 = dataBase(pat).ccep_clin.stimsets_avg(stimp,1);
        stimp2 = dataBase(pat).ccep_clin.stimsets_avg(stimp,2);
          
        % Location of stimulation pair electrodes
        lobe_stimp1 = elec_lobe_label(stimp1,1);
        lobe_stimp2 = elec_lobe_label(stimp2,1);
        
        for m = 1:size(mode,2)
            if isequal(mode{m},'Clin')
                data = dataBase(pat).ccep_clin;
                
            elseif isequal(mode{m},'Prop')
                data = dataBase(pat).ccep_prop;

            end

            % When stimp1 and stimp2 are part of the same lobe
            if isequal(lobe_stimp1{:}, lobe_stimp2{:})

                 if isequal(lobe_stimp1{:},'Temporal') 
                    lat_temporal_stimp.(sprintf('%s',mode{m}))(sz_row, c_temp) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;
                    c_temp = c_temp+1;
    
                elseif isequal(lobe_stimp1{:},'Frontal')
                    lat_frontal_stimp.(sprintf('%s',mode{m}))(sz_row, c_front) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;          
                    c_front = c_front+1;
    
                elseif isequal(lobe_stimp1{:},'Parietal')
                    lat_parietal_stimp.(sprintf('%s',mode{m}))(sz_row, c_pari) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;
                    c_pari = c_pari+1;    
                       
                elseif isequal(lobe_stimp1{:},'Central')
                    lat_central_stimp.(sprintf('%s',mode{m}))(sz_row, c_cent) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;    
                    c_cent = c_cent+1;               
                end

            else % when stimp1 and stimp2 are not on the same lobe
                 
                if isequal(lobe_stimp1{:},'Temporal') || isequal(lobe_stimp2{:},'Temporal')
                    lat_temporal_stimp.(sprintf('%s',mode{m}))(sz_row, c_temp) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;
                    c_temp = c_temp+1;
                end
    
                if isequal(lobe_stimp1{:},'Frontal') || isequal(lobe_stimp2{:},'Frontal')
                    lat_frontal_stimp.(sprintf('%s',mode{m}))(sz_row, c_front) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;          
                    c_front = c_front+1;
    
                end

                if isequal(lobe_stimp1{:},'Parietal') || isequal(lobe_stimp2{:},'Parietal')
                    lat_parietal_stimp.(sprintf('%s',mode{m}))(sz_row, c_pari) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;
                    c_pari = c_pari+1;    
                end

                if isequal(lobe_stimp1{:},'Central') || isequal(lobe_stimp2{:},'Central')
                    lat_central_stimp.(sprintf('%s',mode{m}))(sz_row, c_cent) = (data.n1_peak_sample(:,stimp) -(2*fs)) *ts *1000;    
                    c_cent = c_cent+1;               
                end



            end
            

        end

    end

end 

% Replace zeros with NaN to avoid influence on median
lat_temporal_stimp.Clin(ismember(lat_temporal_stimp.Clin,0)) = NaN;
lat_frontal_stimp.Clin(ismember(lat_frontal_stimp.Clin,0)) = NaN;
lat_parietal_stimp.Clin(ismember(lat_parietal_stimp.Clin,0)) = NaN;
lat_central_stimp.Clin(ismember(lat_central_stimp.Clin,0)) = NaN;

lat_temporal_stimp.Prop(ismember(lat_temporal_stimp.Prop,0)) = NaN;
lat_frontal_stimp.Prop(ismember(lat_frontal_stimp.Prop,0)) = NaN;
lat_parietal_stimp.Prop(ismember(lat_parietal_stimp.Prop,0)) = NaN;
lat_central_stimp.Prop(ismember(lat_central_stimp.Prop,0)) = NaN;

% Make array of latencies in matrix
lat_temporal_stimp.Clin(isnan(lat_temporal_stimp.Clin)) = [];
lat_temporal_stimp.Prop(isnan(lat_temporal_stimp.Prop)) = [];
lat_frontal_stimp.Clin(isnan(lat_frontal_stimp.Clin)) = [];
lat_frontal_stimp.Prop(isnan(lat_frontal_stimp.Prop)) = [];
lat_parietal_stimp.Clin(isnan(lat_parietal_stimp.Clin)) = [];
lat_parietal_stimp.Prop(isnan(lat_parietal_stimp.Prop)) = [];
lat_central_stimp.Clin(isnan(lat_central_stimp.Clin)) = [];
lat_central_stimp.Prop(isnan(lat_central_stimp.Prop)) = [];


%% %% Display all latencies per lobe in a violinplot
violin_mat = NaN(2000,8);   % MAKE THE SIZE OF THIS MATRIX LARGE ENOUGH BECAUSE LARGER ARRAYS ARE FILLED WITH 0 AND INFLUENCE THE MEDIAN
violin_mat(1:size(lat_temporal_stimp.Clin,2), 1) = lat_temporal_stimp.Clin';    % Temporal Clin
violin_mat(1:size(lat_temporal_stimp.Prop,2), 2) = lat_temporal_stimp.Prop';    % Temporal Prop

violin_mat(1:size(lat_frontal_stimp.Clin,2), 3) = lat_frontal_stimp.Clin';    % Frontal Clin
violin_mat(1:size(lat_frontal_stimp.Prop,2), 4) = lat_frontal_stimp.Prop';    % Frontal Prop

violin_mat(1:size(lat_parietal_stimp.Clin,2), 5) = lat_parietal_stimp.Clin';    % Parietal Clin
violin_mat(1:size(lat_parietal_stimp.Prop,2), 6) = lat_parietal_stimp.Prop';    % Parietal Prop

violin_mat(1:size(lat_central_stimp.Clin,2), 7) = lat_central_stimp.Clin';    % Central Clin
violin_mat(1:size(lat_central_stimp.Prop,2), 8) = lat_central_stimp.Prop';    % Central Prop


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


title(sprintf('N1-Latency evoked by stimulation pairs per lobe '),'FontSize', 15, 'FontWeight', 'bold')
ylabel('N1-latency of response evoked per lobe','FontSize', 15, 'FontWeight', 'bold')


% STATISTICS
p_temp = ranksum(lat_temporal_stimp.Clin, lat_temporal_stimp.Prop);
p_fron = ranksum(lat_frontal_stimp.Clin, lat_frontal_stimp.Prop);
p_pari = ranksum(lat_parietal_stimp.Clin, lat_parietal_stimp.Prop);
p_cent = ranksum(lat_central_stimp.Clin, lat_central_stimp.Prop);


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
outlabel=sprintf('N1-latency_lobe_Stimp.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')





end
