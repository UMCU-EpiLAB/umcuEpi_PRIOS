function distance_elec_stimp(dataBase, myDataPath, av_lat_elec)
% This function is used to determine the distance between electrodes and
% the stimulation pair.
% Further on, we can correlate the distance with the N1-peak latency.

% load MNI electrode positions 
load(fullfile(myDataPath.CCEPpath,'elec_coordinatesMNI305.mat'),'elec_coords')

% Preallocate (900 is a raw guess for the total numer of ERs found in a patient.  
distance = NaN(900,size(dataBase,2));
% group = cell(900, size(dataBase,2));
% rcl =1;
% rpr = 1;

% mode = {'clin','prop'};
r = 1;  % used to save the distance-group for all subjects in 1 array

for pat = 1:size(dataBase,2) 

        data = dataBase(pat).ccep_clin; % Does not matter whether you take clin or prop since both sessions should contain the same electrodes
        r_dis = 1;  % save the distance per subject

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
                     % electrodes with a SPES response
                     for elec = 1:size(data.ch,1)
                         % only determine the distance for electrodes-stimpair
                         % combinations that showed an ER. Also exlude
                         % electrodes that are in depth
    
                         if ~isnan(dataBase(pat).ccep_clin.n1_peak_sample(elec, stimp)) && ~isnan(dataBase(pat).ccep_prop.n1_peak_sample(elec, stimp)) && ~isequal(data.tb_channels.group{elec},'depth')
                             % Find name of elec
                             elec_name = data.ch(elec);
                             % Find row number of elec in MNI coordinates matrix
                             loc_elec = find(ismember(table2cell(elec_coords(pat).elecs_tsv(:,1)), elec_name), 1);
            
                             % Coordinates of elec
                             coor_elec = elec_coords(pat).mni_coords(loc_elec,1:3);
            
                             % Distance of elec to middle of the stimpair
                             distance(r_dis,pat) = norm(coor_elec-coor_stimp);        
                             
                             if distance(r_dis,pat) < 20.0
                                 elec_group(r,:) = {'A'};
                             elseif distance(r_dis,pat) > 20.0 && distance(r_dis,pat) < 40.0
                                 elec_group(r,:) = {'B'};
                             elseif distance(r_dis,pat) > 40.0
                                 elec_group(r,:) = {'C'};
                             else
                                 warning(sprintf('went wrong on elec = %d and stimp = %d for pat = %d',elec, stimp, pat))
                             end
                             
                             r_dis = r_dis+1;
                             r = r+1;
            
                         else
                             % do nothing because electrode-stimpair combination
                             % did not result in an ER during clinical-SPES AND propofol-SPES.
                             % OR was a DEPTH electrode
    
                         end        
                     end

            end
        end

end
 
% Remove all nan from distance matrix             

% % Count number of groupmembers
% a=unique(group,'stable');
% b=cellfun(@(x) sum(ismember(group,x)),a,'un',0)

%% Statistics (Wilcoxon rank test)
% Compare medians 
latencyCl = vertcat(dataBase.lat_elec_clin);
latencyPr = vertcat(dataBase.lat_elec_prop);

mode = {'<2','2-4','>4'};
for m = 1:size(mode,2)
    if isequal(mode{m}, '<2')
        group = {'A'};
    elseif isequal(mode{m}, '2-4')
        group = {'B'};
    
    elseif isequal(mode{m}, '>4')
        group = {'C'};
    end

    loc_group = ismember(elec_group, group);
    lat_clin = latencyCl(loc_group);
    lat_prop = latencyPr(loc_group);
    
    p_group = signrank(lat_clin, lat_prop);

    fprintf('P-value between clin and prop for group with %s cm = %1.3f \n', mode{m}, p_group)
    fprintf('Median group with %s, clin = %1.2f, prop = %1.2f \n', mode{m}, median(lat_clin), median(lat_prop))

end


%% Compare median latency of electrodes in SOZ and not-SOZ between propofol and clin

% 100 rows is a rough estimation! make sure to check this otherwise it is 
% filled WITH ZEROS and these are not ignored in the calculations
lat_soz_clin = NaN(100,size(dataBase,2));  
lat_non_soz_clin = NaN(100,size(dataBase,2));
lat_soz_prop = NaN(100,size(dataBase,2));
lat_non_soz_prop = NaN(100,size(dataBase,2));

for pat = 1:size(elec_coords,2)

    soz_elec = ismember(elec_coords(pat).elecs_tsv.soz  ,'yes');
    soz_name = elec_coords(pat).elecs_tsv(soz_elec,1);
    soz_elec_loc = ismember(dataBase(pat).ccep_clin.ch, table2cell(soz_name));   % does not matter whether you take clin or prop, channels are equal       
    non_soz_elec_loc = ~ismember(dataBase(pat).ccep_clin.ch, table2cell(soz_name));
    
    if sum(soz_elec) > 1        % only plot violin when more than 1 
        % av_lat_elec (average latency per electrode) is already calculated in
        % function plot_electrodes_on_MRI
        clin_colm = pat*2-1;
        prop_colm = pat*2;
    
        lat_soz_clin(1:sum(soz_elec_loc),pat) = av_lat_elec(soz_elec_loc, clin_colm);
        lat_non_soz_clin(1:sum(non_soz_elec_loc),pat) = av_lat_elec(non_soz_elec_loc, clin_colm);
    
        lat_soz_prop(1:sum(soz_elec_loc),pat) = av_lat_elec(soz_elec_loc, prop_colm);
        lat_non_soz_prop(1:sum(non_soz_elec_loc),pat) = av_lat_elec(non_soz_elec_loc, prop_colm);
    
        % When variables in the set are different (non-paired) than use rank
        % sum 
        % --> soz with non-soz
        p_clin = ranksum(lat_soz_clin(:,pat), lat_non_soz_clin(:,pat));
        p_prop = ranksum(lat_soz_prop(:,pat), lat_non_soz_prop(:,pat));
       
        % When comparing two measurements of the same data set than use
        % signrank. 
        % --> soz_clin with soz_prop
        p_soz = signrank(lat_soz_clin(:,pat), lat_soz_prop(:,pat));
        p_non_soz = signrank(lat_non_soz_clin(:,pat), lat_non_soz_prop(:,pat));
        
        %preallocatio 
        violin_mat = NaN(size(lat_non_soz_clin(:,pat),1),4);
        violin_mat(1:size(lat_soz_clin(:,pat),1),1) = lat_soz_clin(:,pat);
        violin_mat(1:size(lat_non_soz_clin(:,pat),1),2) = lat_non_soz_clin(:,pat);
        violin_mat(1:size(lat_soz_prop(:,pat),1),3) = lat_soz_prop(:,pat);
        violin_mat(1:size(lat_non_soz_prop(:,pat),1),4) = lat_non_soz_prop(:,pat);    
       
        figure('Position',[205,424,1530,638]);
        violins = violinplot(violin_mat) ;                    
        for i = 1:2:size(violin_mat,2)
            violins(i).ViolinColor(:) = [1 0 0];
            violins(i+1).ViolinColor(:) = [0 0 1];
        end
    
        ax = gca;
        ax.XTickLabel = [{'clinical SOZ'} {'clinical non-SOZ'} {'propofol SOZ'} {'propofol non-SOZ'}];
    
        title(sprintf('N1-latency of electrodes in or outside SOZ for subj %s', dataBase(pat).ccep_clin.sub_label),'FontSize', 15, 'FontWeight', 'bold')
        ylabel('N1-latency per electrode','FontSize', 15, 'FontWeight', 'bold')
    
        medians = median(violin_mat,'omitnan');
        ymin = min(ylim);
        y_range = diff(ylim);
        x_as = 1:size(violin_mat,2);
        size_pat = size(violin_mat,2); 
        second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
        text([(x_as(1)-x_as(2))*-0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')
    
    
         % Set significance in plot with a line and not in name
         count = 1;
         ymax = max(max(violin_mat));
         
         if p_clin  < 0.01 
             text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
             plot(count+0.01:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    
         elseif p_clin  < 0.05 
             text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
         end
    
         if p_soz  < 0.01 
             text(count+1,ymax-5,'**','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)
    
    
         elseif p_soz  < 0.05 
             text(count+1,ymax-5,'*','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+1.9, ymax-5.93*ones(19,1),'k','LineWidth',2)
    
          end
    
    
        count = 3;
         if p_prop  < 0.01 
             text(count+0.5,ymax-10,'**','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
    
         elseif p_prop  < 0.05 
             text(count+0.5,ymax-10,'*','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+0.9, ymax-10.93*ones(9,1),'k','LineWidth',2)
    
         end
    
         count =2 ;
          if p_non_soz  < 0.01 
             text(count+1,ymax-0.1,'**','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)
    
    
         elseif p_non_soz  < 0.05 
             text(count+1,ymax-0.1,'*','FontSize',20,'FontWeight','bold')
             plot(count+0.1:0.1:count+1.9, ymax-0.93*ones(19,1),'k','LineWidth',2)
    
          end
    end

end

%% plot for all patients combined
    % Determine the correlation of latency and soz or non-soz when grouping all patients
    p_clin_all = ranksum(lat_soz_clin(:), lat_non_soz_clin(:));
    p_prop_all = ranksum(lat_soz_prop(:), lat_non_soz_prop(:));
    
    p_soz_all = signrank(lat_soz_clin(:), lat_soz_prop(:));
    p_non_soz_all = signrank(lat_non_soz_clin(:), lat_non_soz_prop(:));
    
    fprintf('P-value between SOZ and non-SOZ for for all patients during clinical = %1.3f \n',p_clin_all);
    fprintf('P-value between SOZ and non-SOZ for for all patients during propofol = %1.3f \n',p_prop_all);
    fprintf('P-value between SOZ-clin and SOZ-prop for for all patients = %1.3f \n',p_soz_all);
    fprintf('P-value between non-SOZ-clin and non-SOZ-prop for for all patients = %1.3f \n',p_non_soz_all);
    
    %preallocatio 
    violin_mat = NaN(size(lat_non_soz_clin(:),1),4);
    violin_mat(:,1) = lat_soz_clin(:);
    violin_mat(:,2) = lat_non_soz_clin(:);
    violin_mat(:,3) = lat_soz_prop(:);
    violin_mat(:,4) = lat_non_soz_prop(:);    
    
    figure('Position',[205,424,1530,638]);
    violins = violinplot(violin_mat) ;                    
    violins(1).ViolinColor(:) = [1 0 0]; violins(2).ViolinColor(:) = [1 0 0];
    violins(3).ViolinColor(:) = [0 0 1]; violins(4).ViolinColor(:) = [0 0 1];    
    
    ax = gca;
    ax.XTickLabel = [{'clinical SOZ'} {'clinical non-SOZ'} {'propofol SOZ'} {'propofol non-SOZ'}];
    
    title(sprintf('N1-latency of electrodes in or outside SOZ for all patients'),'FontSize', 15, 'FontWeight', 'bold')
    ylabel('N1-latency per electrode','FontSize', 15, 'FontWeight', 'bold')
    
    medians = median(violin_mat,'omitnan');
    ymin = min(ylim);
    y_range = diff(ylim);
    x_as = 1:size(violin_mat,2);
    size_pat = size(violin_mat,2); 
    second_row_txt = cellstr(strsplit(num2str(medians,'%.2f '),' '));
    text([(x_as(1)-x_as(2))*-0.5 x_as], ones(1,size_pat+1)*ymin-0.1*y_range, ['Median' second_row_txt],'HorizontalAlignment','center','FontSize', 12, 'FontWeight', 'bold')
    
    % Set significance in plot with a line and not in name
    count = 1;
    ymax = max(max(violin_mat));
    
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
outlabel=sprintf('N1-latency_soz_non-SOZ.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')


end