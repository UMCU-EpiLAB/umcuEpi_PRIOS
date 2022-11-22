function distance_elec_stimp(dataBase, myDataPath)
% This function is used to determine the distance between electrodes and
% the stimulation pair.
% Further on, we can correlate the distance with the N1-peak latency.

% load MNI electrode positions 
load(fullfile(myDataPath.CCEPpath,'elec_coordinatesMNI305.mat'),'elec_coords')

fs = 2048;
ts = 1/fs ;

r = 1;  % used to save the distance-group for all subjects in 1 array
r_clin = 1;
r_prop = 1;

for pat = 1:size(dataBase,2) 

    data = dataBase(pat).ccep_clin; % Does not matter whether you take clin or prop since both sessions should contain the same electrodes
%         r_dis = 1;  % save the distance per subject

    for stimp = 1:size(data.stimpnames_avg,2)
        % Determine the two electrodes in the stimulation pair
        elec1 = extractBefore(data.stimpnames_avg(1,stimp),'-');
        elec2 = extractAfter(data.stimpnames_avg(1,stimp),'-');

        % Find the location of the electrodes in the matrix
        loc1 = find(ismember(table2cell(elec_coords(pat).elecs_tsv(:,1)), elec1), 1);
        loc2 = find(ismember(table2cell(elec_coords(pat).elecs_tsv(:,1)), elec2), 1);

        if isempty(loc1) || isempty(loc2)
            error('Electrode from stimulation pair is not found in coordinates matrix')
        end

        % If electrode is an depth electrode, exclude from the analysis
        % since the distance between depth electrodes is different from
        % the distance between coretx electrodes
        if ismember(elec_coords(pat).elecs_tsv.group(loc1),'depth') || ismember(elec_coords(pat).elecs_tsv.group(loc2),'depth')
            % Skip during the distance analysis
        else
            % Find the coordinates of the stimpair electrodes
            coor1 = elec_coords(pat).mni_coords(loc1,1:3);
            coor2 = elec_coords(pat).mni_coords(loc2,1:3);
    
            % Find the distance between the two electrodes of the stimpair
            coor_stimp(:,1) = (coor1(1)+coor2(1))/2;
            coor_stimp(:,2) = (coor1(2)+coor2(2))/2;
            coor_stimp(:,3) = (coor1(3)+coor2(3))/2;
       
    
             % Check that the middle between the two stimpair electrodes is correct
             if ~isequal(norm(coor1-coor_stimp), norm(coor_stimp-coor2))
                 error('did not find the correct middle of the stimulation pair electrodes, row 59')
             end
                        
                 % Determine distance between stimpair coor and all other
                 % electrodes with a SPES respons. EXCLUDE DEPTH
                 % ELECTRODES   
                 for elec = 1:size(data.ch,1)

                     % Find name of elec
                     elec_name = data.ch(elec);
                     % Find row number of elec in MNI coordinates matrix
                     loc_elec = find(ismember(table2cell(elec_coords(pat).elecs_tsv(:,1)), elec_name), 1);
    
                     % Coordinates of elec
                     coor_elec = elec_coords(pat).mni_coords(loc_elec,1:3);
                      
                     % Distance of elec to middle of the stimpair
                     distance = norm(coor_elec-coor_stimp); 

%                          if ~isnan(dataBase(pat).ccep_clin.n1_peak_sample(elec, stimp)) && ~isnan(dataBase(pat).ccep_prop.n1_peak_sample(elec, stimp)) && ~isequal(data.tb_channels.group{elec},'depth')                                                                         
%                              
%                              if distance < 20.0
%                                  elec_group(r,:) = {'A'};
%                              elseif distance > 20.0 && distance < 40.0
%                                  elec_group(r,:) = {'B'};
%                              elseif distance > 40.0
%                                  elec_group(r,:) = {'C'};
%                              else
%                                  warning(sprintf('Something went wrong on elec = %d and stimp = %d for pat = %d',elec, stimp, pat))
%                              end
%                              
%                              r = r+1;
%             
%                          else
%                              % do nothing because electrode-stimpair combination
%                              % did not result in an ER during clinical-SPES AND propofol-SPES.
%                              % OR was a DEPTH electrode
%     
%                          end 


                     %% For responses in clinical-SPES
                     % Independent of the response in propofol-SPES
                     if ~isnan(dataBase(pat).ccep_clin.n1_peak_sample(elec, stimp)) && ~isequal(data.tb_channels.group{elec},'depth')

%                              if distance < 20.0
%                                  dist_clin(r_clin,pat) = {'A'};
%                              elseif distance > 20.0 && distance < 40.0
%                                  dist_clin(r_clin,pat) = {'B'};
%                              elseif distance > 40.0
%                                  dist_clin(r_clin,pat) = {'C'};
%                              else
%                                  warning(sprintf('SOmething went wrong on elec = %d and stimp = %d for pat = %d',elec, stimp, pat))
%                              end
                         

                         % Save the absolute number of distance
                         % between stimpair and electrode and the
                         % latency in the second column
                         dist_lat_clin(r_clin,1) = distance;
                         lat_clin = dataBase(pat).ccep_clin.n1_peak_sample(elec, stimp) - (2*fs);
                         dist_lat_clin(r_clin,2) = lat_clin *ts *1000;    % in ms

                         r_clin = r_clin+1;                             

                     end

                     % For responses in propofol-SPES
                     if ~isnan(dataBase(pat).ccep_prop.n1_peak_sample(elec, stimp)) && ~isequal(data.tb_channels.group{elec},'depth')
                                                      
%                              if distance < 20.0
%                                  dist_prop(r_prop,:) = {'A'};
%                              elseif distance > 20.0 && distance < 40.0
%                                  dist_prop(r_prop,:) = {'B'};
%                              elseif distance > 40.0
%                                  dist_prop(r_prop,:) = {'C'};
%                              else
%                                  warning(sprintf('went wrong on elec = %d and stimp = %d for pat = %d',elec, stimp, pat))
%                              end


                         % Save the absolute number of distance
                         % between stimpair and electrode and the
                         % latency in the second column
                         dist_lat_prop(r_prop,1) = distance;
                         lat_prop = dataBase(pat).ccep_prop.n1_peak_sample(elec, stimp) - (2*fs);
                         dist_lat_prop(r_prop,2) = lat_prop *ts *1000;    % in ms
                         
                         r_prop = r_prop +1;                                                   

                     end 

                 end

        end
    end

end



% Make curve with distance and latency
% Use circle for clinical-SPES results and squares for propofol-SPES
figure('Position',[393,540,847,513])
coefficients_clin = polyfit(dist_lat_clin(:,1),dist_lat_clin(:,2), 1);
xFit_clin = linspace(min(dist_lat_clin(:,1)), max(dist_lat_clin(:,1)), 2000); % Option 2 : lots of points, and not just where the training points are.
yFit_clin = polyval(coefficients_clin, xFit_clin);
scatter(dist_lat_clin(:,1), dist_lat_clin(:,2), 'o', 'MarkerEdgeColor','r','MarkerEdgeAlpha',0.2)
hold on;

coefficients = polyfit(dist_lat_prop(:,1),dist_lat_prop(:,2), 1);
xFit = linspace(min(dist_lat_prop(:,1)), max(dist_lat_prop(:,1)), 2000); % Option 2 : lots of points, and not just where the training points are.
yFit = polyval(coefficients, xFit);
scatter(dist_lat_prop(:,1), dist_lat_prop(:,2), '*','MarkerEdgeColor','b','MarkerEdgeAlpha',0.2)
hold on
plot(xFit_clin, yFit_clin, 'LineWidth', 3,'Color',[0.64,0.08,0.18])
plot(xFit, yFit, 'LineWidth', 3,'Color',[0.00,0.45,0.74])


legend('Clinical-SPES','Propofol-SPES','Fit for Clinical-SPES','Fit for Propofol-SPES')

xlabel('Distance between stimulation pair and electrode (mm)')
ylabel('Latency of N1-peak (ms)')


% Statistics (wilcoxon rank sum test)
% non-paired, not normally distributed (non-parametric)
Latency = [dist_lat_clin(:,2); dist_lat_prop(:,2)];
Distance = [dist_lat_clin(:,1); dist_lat_prop(:,1)];

Session = cell(size(Distance,1),1);
Session(1:size(dist_lat_clin(:,1),1),:) = {'clin'};
Session(size(dist_lat_clin(:,1),1)+1:end,:) = {'prop'};

% When ANOVA table shows that Distance and Session correlation is <0.05
% there is a significant differnce betweeen the slopes caused by the
% distance

% When line is steeper --> The independent variable has a stronger effect on the dependent variable
% in one categoriythan it does to the other category
[h,atab,ctab,stats] = aoctool(Distance,Latency,Session,0.05);

%%
 
% Remove all nan from distance matrix             

% % Count number of groupmembers
% a=unique(group,'stable');
% b=cellfun(@(x) sum(ismember(group,x)),a,'un',0)


%% Create array with the number of ERs evoked per electrode for ALL patients
for pat = 1:size(dataBase,2)
    dataBase(pat).ERs_per_elec_clin = dataBase(pat).agreement_parameter.ERs_elecClin';
    dataBase(pat).ERs_per_elec_prop = dataBase(pat).agreement_parameter.ERs_elecProp';
end

ERs_per_elec_clin = vertcat(dataBase.ERs_per_elec_clin);
ERs_per_elec_prop = vertcat(dataBase.ERs_per_elec_prop);


%% Statistics (Wilcoxon rank test)
% Compare medians 
% latencyCl = vertcat(dataBase.lat_elec_clin);
% latencyPr = vertcat(dataBase.lat_elec_prop);
% 
% mode = {'<2','2-4','>4'};
% for m = 1:size(mode,2)
%     if isequal(mode{m}, '<2')
%         group = {'A'};
%     elseif isequal(mode{m}, '2-4')
%         group = {'B'};
%     
%     elseif isequal(mode{m}, '>4')
%         group = {'C'};
%     end
% 
%     loc_group = ismember(elec_group, group);
%     lat_clin = latencyCl(loc_group);
%     lat_prop = latencyPr(loc_group);
%     
%     p_group = signrank(lat_clin, lat_prop);
% 
%     fprintf('P-value between clin and prop for group with %s cm = %1.3f \n', mode{m}, p_group)
%     fprintf('Median group with %s, clin = %1.2f, prop = %1.2f \n', mode{m}, median(lat_clin), median(lat_prop))
% 
% end


%% Compare median latency of electrodes in SOZ and not-SOZ between propofol and clin
% 
% % 100 rows is a rough estimation! make sure to check this otherwise it is 
% % filled WITH ZEROS and these are not ignored in the calculations
% lat_soz_clin = NaN(100,size(dataBase,2));  
% lat_non_soz_clin = NaN(100,size(dataBase,2));
% lat_soz_prop = NaN(100,size(dataBase,2));
% lat_non_soz_prop = NaN(100,size(dataBase,2));

for pat = 1:size(elec_coords,2)

    soz_elec = ismember(elec_coords(pat).elecs_tsv.soz  ,'yes');
    soz_name = elec_coords(pat).elecs_tsv(soz_elec,1);
    soz_elec_loc = ismember(dataBase(pat).ccep_clin.ch, table2cell(soz_name));   % does not matter whether you take clin or prop, channels are equal       
%     non_soz_elec_loc = ~ismember(dataBase(pat).ccep_clin.ch, table2cell(soz_name));
    
%     if sum(soz_elec) > 1        % only plot violin when more than 1 
       
    % Calculate median N1-latency of all responses evoked on electrode on
    % the SOZ and outside
    idx_elec = find(soz_elec_loc == 1);
    idx_nelec = find(soz_elec_loc == 0);
   
    if pat == 1
        row_start = 1;
        row_start_nS = 1;
    else 
        row_start = size(lat_elec_SOZ,1)+1;
        row_start_nS = size(lat_elec_nSOZ,1)+1;
    end

    lat_elec_SOZ(row_start:(row_start+size(idx_elec,1)-1), 1:size(dataBase(pat).ccep_clin.n1_peak_sample,2)) = dataBase(pat).ccep_clin.n1_peak_sample(idx_elec,:);
    lat_elec_nSOZ(row_start_nS:(row_start_nS+size(idx_nelec,1)-1), 1:size(dataBase(pat).ccep_clin.n1_peak_sample,2)) = dataBase(pat).ccep_clin.n1_peak_sample(idx_nelec,:);
    
    lat_elec_SOZ_prop(row_start:(row_start+size(idx_elec,1)-1), 1:size(dataBase(pat).ccep_prop.n1_peak_sample,2)) = dataBase(pat).ccep_prop.n1_peak_sample(idx_elec,:);
    lat_elec_nSOZ_prop(row_start_nS:(row_start_nS+size(idx_nelec,1)-1), 1:size(dataBase(pat).ccep_prop.n1_peak_sample,2)) = dataBase(pat).ccep_prop.n1_peak_sample(idx_nelec,:);
         
        
      
%%  Determine if propofol has effect on the N1-latency when one of the stimulation pair electrodes is located on the SOZ
    Stimp_on_SOZ = zeros(size(dataBase(pat).ccep_clin.stimpnames_avg,2),1);

    for stimp = 1:size(dataBase(pat).ccep_clin.stimpnames_avg,2)
        % Split names to obtain electrodes in stimulation pair
        stimpname1 = extractBefore(dataBase(pat).ccep_clin.stimpnames_avg(stimp),'-');
        stimpname2 = extractAfter(dataBase(pat).ccep_clin.stimpnames_avg(stimp),'-');

        % If either one of the stimulation pair electrodes is on the SOZ,
        % than place it in the SOZ group.
        if ismember(stimpname1, table2cell(soz_name)) || ismember(stimpname2, table2cell(soz_name))
            Stimp_on_SOZ(stimp) = 1;
        else
            Stimp_on_SOZ(stimp) = 0;
        end        

    end

    % Calculate median N1-latency of all responses evoked by stimulation
    % pairs on the SOZ
    idx_stimp = find(Stimp_on_SOZ == 1);
    idx_nstimp = find(Stimp_on_SOZ == 0);
   
    if pat == 1
        col_start = 1;
        col_start_nS = 1;
    else 
        col_start = size(lat_stimp_SOZ,2)+1;
        col_start_nS = size(lat_stimp_nSOZ,2)+1;
    end

    lat_stimp_SOZ(1:size(dataBase(pat).ccep_clin.n1_peak_sample,1),col_start:(col_start+size(idx_stimp,1)-1)) = dataBase(pat).ccep_clin.n1_peak_sample(:,idx_stimp);
    lat_stimp_nSOZ(1:size(dataBase(pat).ccep_clin.n1_peak_sample,1),col_start_nS:(col_start_nS+size(idx_nstimp,1)-1)) = dataBase(pat).ccep_clin.n1_peak_sample(:,idx_nstimp);
    
    lat_stimp_SOZ_prop(1:size(dataBase(pat).ccep_prop.n1_peak_sample,1),col_start:(col_start+size(idx_stimp,1)-1)) = dataBase(pat).ccep_prop.n1_peak_sample(:,idx_stimp);
    lat_stimp_nSOZ_prop(1:size(dataBase(pat).ccep_prop.n1_peak_sample,1),col_start_nS:(col_start_nS+size(idx_nstimp,1)-1)) = dataBase(pat).ccep_prop.n1_peak_sample(:,idx_nstimp);
    
    
end

% Convert zeros to NaN otherwise they influence the median
lat_elec_SOZ(ismember(lat_elec_SOZ,0)) = NaN;
lat_elec_nSOZ(ismember(lat_elec_nSOZ,0)) = NaN;    
lat_elec_SOZ_prop(ismember(lat_elec_SOZ_prop,0)) = NaN;
lat_elec_nSOZ_prop(ismember(lat_elec_nSOZ_prop,0)) = NaN;

lat_stimp_SOZ(ismember(lat_stimp_SOZ,0)) = NaN;
lat_stimp_nSOZ(ismember(lat_stimp_nSOZ,0)) = NaN;
lat_stimp_SOZ_prop(ismember(lat_stimp_SOZ_prop,0)) = NaN;
lat_stimp_nSOZ_prop(ismember(lat_stimp_nSOZ_prop,0)) = NaN;

% Calculate the median of the N1-latency of the SOZ stimulation pairs and
% the Non-SOZ stimulation pairs
lat_elec_SOZ_sec = (lat_elec_SOZ-(2*fs))*ts*1000;  
lat_elec_nSOZ_sec = (lat_elec_nSOZ-(2*fs))*ts*1000;  
lat_elec_SOZ_prop_sec = (lat_elec_SOZ_prop-(2*fs))*ts*1000;  
lat_elec_nSOZ_prop_sec = (lat_elec_nSOZ_prop-(2*fs))*ts*1000;  


lat_SOZ_stimp_sec = (lat_stimp_SOZ-(2*fs))*ts*1000;  
lat_nSOZ_stimp_sec = (lat_stimp_nSOZ-(2*fs))*ts*1000;  
lat_SOZ_stimp_sec_prop = (lat_stimp_SOZ_prop-(2*fs))*ts*1000;  
lat_nSOZ_stimp_sec_prop = (lat_stimp_nSOZ_prop-(2*fs))*ts*1000;  


%% Display all latencies of SOZ and non-SOZ electrodes in a violinplot
violin_mat = NaN(64000,4);
violin_mat(1:size(lat_elec_SOZ_sec,1)*size(lat_elec_SOZ_sec,2),1) = vertcat(lat_elec_SOZ_sec(:));
violin_mat(1:size(lat_elec_nSOZ_sec,1)*size(lat_elec_nSOZ_sec,2),2) = vertcat(lat_elec_nSOZ_sec(:));
violin_mat(1:size(lat_elec_SOZ_prop_sec,1)*size(lat_elec_SOZ_prop_sec,2),3) = vertcat(lat_elec_SOZ_prop_sec(:));
violin_mat(1:size(lat_elec_nSOZ_prop_sec,1)*size(lat_elec_nSOZ_prop_sec,2),4) = vertcat(lat_elec_nSOZ_prop_sec(:));    

figure('Position',[205,424,1530,638]);
violins = violinplot(violin_mat) ;                    
violins(1).ViolinColor(:) = [1 0 0]; violins(2).ViolinColor(:) = [1 0 0];
violins(3).ViolinColor(:) = [0 0 1]; violins(4).ViolinColor(:) = [0 0 1];  

% Calculate the number of CCEPs evoked during each session
nr_soz = sum(~isnan(violin_mat(:,1)));
nr_nsoz = sum(~isnan(violin_mat(:,2)));
nr_soz_p = sum(~isnan(violin_mat(:,3)));
nr_nsoz_p = sum(~isnan(violin_mat(:,4)));

ax = gca;
ax.XTickLabel = [{sprintf('Clinical on SOZ, n = %d',nr_soz)}  {sprintf('Clinical not on SOZ, n = %d', nr_nsoz)} ...
            {sprintf('Propofol on SOZ, n = %d',nr_soz_p)} {sprintf('Propofol not on SOZ, n = %d', nr_nsoz_p)}];

title(sprintf('N1-latency of electrodes on or outside SOZ for all patients'),'FontSize', 15, 'FontWeight', 'bold')
ylabel('N1-latency of response detected per electrode','FontSize', 15, 'FontWeight', 'bold')

medians = median(violin_mat,'omitnan');
ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(violin_mat,2);
size_pat = size(violin_mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
text([(x_as(1)-x_as(2))*-0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')

% Make array of latencies in matrix and remove NaNs
lat_elec_SOZ_sec(isnan(lat_elec_SOZ_sec))= [];
lat_elec_nSOZ_sec(isnan(lat_elec_nSOZ_sec)) = [];
lat_elec_SOZ_prop_sec(isnan(lat_elec_SOZ_prop_sec)) = [];
lat_elec_nSOZ_prop_sec(isnan(lat_elec_nSOZ_prop_sec)) = [];

% STATISTICS
p_clin_all = ranksum(lat_elec_SOZ_sec(:), lat_elec_nSOZ_sec(:));
p_prop_all = ranksum(lat_elec_SOZ_prop_sec(:), lat_elec_nSOZ_prop_sec(:));

p_soz_all = ranksum(lat_elec_SOZ_sec(:), lat_elec_SOZ_prop_sec(:));
p_non_soz_all = ranksum(lat_elec_nSOZ_sec(:), lat_elec_nSOZ_prop_sec(:));

fprintf('ELEC: P-value between SOZ and non-SOZ for for all patients during clinical = %1.3f \n',p_clin_all);
fprintf('ELEC: P-value between SOZ and non-SOZ for for all patients during propofol = %1.3f \n',p_prop_all);
fprintf('ELEC: P-value between SOZ-clin and SOZ-prop for for all patients = %1.3f \n',p_soz_all);
fprintf('ELEC: P-value between non-SOZ-clin and non-SOZ-prop for for all patients = %1.3f \n',p_non_soz_all);


    ymax = max(max(violin_mat));
    count = 1;
    if p_clin_all  < 0.01 
     text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
     plot(count+0.01:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    
    elseif p_clin_all  < 0.05 
     text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    end
    
    if p_soz_all  < 0.01 
     text(count+1,ymax-5,'**','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)
    
    
    elseif p_soz_all  < 0.05 
     text(count+1,ymax-5,'*','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)
    
    end
    
    
    count = 3;
    if p_prop_all  < 0.01 
     text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    
    elseif p_prop_all  < 0.05 
     text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    end
    
    count =2 ;
    if p_non_soz_all  < 0.01 
     text(count+1,ymax-0.1,'**','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)
    
    
    elseif p_non_soz_all  < 0.05 
     text(count+1,ymax-0.1,'*','FontSize',20,'FontWeight','bold')
     plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)
    
    end

    
% Save figure
outlabel=sprintf('ELEC-N1-latency_soz_non-SOZ.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')



%% Display all latencies of SOZ and non-SOZ stimulation pairs in a violinplot
violin_mat = NaN(64000,4);
violin_mat(1:size(lat_SOZ_stimp_sec,1)*size(lat_SOZ_stimp_sec,2),1) = vertcat(lat_SOZ_stimp_sec(:));
violin_mat(1:size(lat_nSOZ_stimp_sec,1)*size(lat_nSOZ_stimp_sec,2),2) = vertcat(lat_nSOZ_stimp_sec(:));
violin_mat(1:size(lat_SOZ_stimp_sec_prop,1)*size(lat_SOZ_stimp_sec_prop,2),3) = vertcat(lat_SOZ_stimp_sec_prop(:));
violin_mat(1:size(lat_nSOZ_stimp_sec_prop,1)*size(lat_nSOZ_stimp_sec_prop,2),4) = vertcat(lat_nSOZ_stimp_sec_prop(:));    

figure('Position',[205,424,1530,638]);
violins = violinplot(violin_mat) ;                    
violins(1).ViolinColor(:) = [1 0 0]; violins(2).ViolinColor(:) = [1 0 0];
violins(3).ViolinColor(:) = [0 0 1]; violins(4).ViolinColor(:) = [0 0 1];  

% Calculate the number of CCEPs evoked during each session
nr_soz = sum(~isnan(violin_mat(:,1)));
nr_nsoz = sum(~isnan(violin_mat(:,2)));
nr_soz_p = sum(~isnan(violin_mat(:,3)));
nr_nsoz_p = sum(~isnan(violin_mat(:,4)));

ax = gca;
ax.XTickLabel = [{sprintf('Clinical on SOZ, n = %d',nr_soz)}  {sprintf('Clinical not on SOZ, n = %d', nr_nsoz)} ...
            {sprintf('Propofol on SOZ, n = %d',nr_soz_p)} {sprintf('Propofol not on SOZ, n = %d', nr_nsoz_p)}];

title(sprintf('N1-latency of stimulation pairs on or outside SOZ for all patients'),'FontSize', 15, 'FontWeight', 'bold')
ylabel('N1-latency of response evoked per stimulation pair','FontSize', 15, 'FontWeight', 'bold')

medians = median(violin_mat,'omitnan');
ymin = min(ylim);
y_range = diff(ylim);
x_as = 1:size(violin_mat,2);
size_pat = size(violin_mat,2); 
second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
text([(x_as(1)-x_as(2))*-0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')
   
% Make array of latencies in matrix and remove NaNs
lat_SOZ_stimp_sec(isnan(lat_SOZ_stimp_sec))= [];
lat_nSOZ_stimp_sec(isnan(lat_nSOZ_stimp_sec)) = [];
lat_SOZ_stimp_sec_prop(isnan(lat_SOZ_stimp_sec_prop)) = [];
lat_nSOZ_stimp_sec_prop(isnan(lat_nSOZ_stimp_sec_prop)) = [];



%% plot for all patients combined
% Determine the correlation of latency and soz or non-soz when grouping all patients
p_clin_all = ranksum(lat_SOZ_stimp_sec(:), lat_nSOZ_stimp_sec(:));
p_prop_all = ranksum(lat_SOZ_stimp_sec_prop(:), lat_nSOZ_stimp_sec_prop(:));

p_soz_all = ranksum(lat_SOZ_stimp_sec(:), lat_SOZ_stimp_sec_prop(:));
p_non_soz_all = ranksum(lat_nSOZ_stimp_sec(:), lat_nSOZ_stimp_sec_prop(:));

fprintf('STIMP: P-value between SOZ and non-SOZ for for all patients during clinical = %1.3f \n',p_clin_all);
fprintf('STIMP: P-value between SOZ and non-SOZ for for all patients during propofol = %1.3f \n',p_prop_all);
fprintf('STIMP: P-value between SOZ-clin and SOZ-prop for for all patients = %1.3f \n',p_soz_all);
fprintf('STIMP: P-value between non-SOZ-clin and non-SOZ-prop for for all patients = %1.3f \n',p_non_soz_all);

ymax = max(max(violin_mat));
count = 1;
   
if p_clin_all  < 0.01 
 text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.01:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)


elseif p_clin_all  < 0.05 
 text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)

end

if p_soz_all  < 0.01 
 text(count+1,ymax-5,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)


elseif p_soz_all  < 0.05 
 text(count+1,ymax-5,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)

end


count = 3;
if p_prop_all  < 0.01 
 text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)


elseif p_prop_all  < 0.05 
 text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)

end

count =2 ;
if p_non_soz_all  < 0.01 
 text(count+1,ymax-0.1,'**','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)


elseif p_non_soz_all  < 0.05 
 text(count+1,ymax-0.1,'*','FontSize',20,'FontWeight','bold')
 plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)

end


% Save figure
outlabel=sprintf('STIMP-N1-latency_soz_non-SOZ.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')


end