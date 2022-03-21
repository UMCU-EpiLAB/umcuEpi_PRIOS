function latency_N1_MRI(av_lat_elec, Lmnipial_vert, Rmnipial_vert, Lmnipial_face, Rmnipial_face, allmni_coords, all_hemi, dataBase, myDataPath)

%% LEFT
% make table with electrodes and the average latency for that electrode
v_d = [270 0];

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      
new_els = els(~isnan(els(:,1)),:);

all_hemi = all_hemi(~isnan(els(:,1)));
start_row = zeros(size(dataBase,2)+1,1);

mode = {'SPESclin','SPESprop'};

for pat = 1:size(dataBase,2)
    % els is nu nog een lange kolom met alle patienten, die moet gesplitst
    % worden
    if pat == 1
        start_row(pat,:) = 1;
    else
        start_row(pat,:) = size(dataBase(pat-1).agreement_parameter.ERs_elecClin,2) + start_row(pat-1,:);
    end
end
start_row(end,:) = size(new_els,1);

for m = 1:size(mode,2)

    figure
    gl.faces = Lmnipial_face+1;
    gl.vertices = Lmnipial_vert;
    gl = gifti(gl);
    tH = ieeg_RenderGifti(gl); %#ok<NASGU>

    for pat = 1:size(dataBase,2)
        
        clin_colm = 2*pat-1;                      
        prop_colm = 2*pat; 
    
        if isequal(mode{m},'SPESclin')
            ERs_elec = av_lat_elec(:,clin_colm);
        elseif isequal(mode{m},'SPESprop')
            ERs_elec = av_lat_elec(:,prop_colm);
        end
    
        if all(ismember(all_hemi(start_row(pat):start_row(pat+1),:),'L'))      
            % make table with number of N1's and coordinates of all electrodes
            pat_elec_in_els = start_row(pat,:) : start_row(pat,:)+size(dataBase(pat).agreement_parameter.ERs_elecClin,2)-1 ;
            els_with_N1 = [new_els(pat_elec_in_els,:) ,ERs_elec(1:size(dataBase(pat).agreement_parameter.ERs_elecClin,2),:)];
            [~,idx] = sort(els_with_N1(:,4),'descend');       % Rank/sort based on the 4th column, high to low
            els_ranked = els_with_N1(idx,:);                               % Rank coordinates as well based on number of ERs      
                   
           % The darkest color indicates the highest latency for each patient. 
           % Colormap
            cm = colormap(flipud(hot(size(unique(els_ranked(:,4)),1)+1)));
            cb = colorbar;
            tickLabels = round(min(unique(els_ranked(:,4))): max(unique(els_ranked(:,4)))/(size(cb.Ticks,2)-1) :max(unique(els_ranked(:,4))));
            
            set(cb, 'TickLabels', cellstr(num2str((tickLabels'))));
            
            unique_color = 1;
            for elec = 1:size(els_ranked,1)
            
                if elec == 1                        % Only the first elec has to start with a unique color
                    ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                    if els_ranked(elec,4) == els_ranked(elec+1,4)
                        % do nothing, unique_color should remain the same
                    elseif els_ranked(elec,4) ~= els_ranked(elec+1,4)
                        unique_color = unique_color + 1; 
                    end
            
                else
                    if els_ranked(elec,4) == els_ranked(elec-1,4) % if the next stimpair has the same number of N1's as the previous, then give the same color
                        ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                    
                    else
                        ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                        unique_color = unique_color + 1;        % if the the number of N1's is different from the previous elec, then go to next color
                    end
                end
    
                hold on
            
            end
        
        else 
            % Do nothing because electrodes of this patient are on the right
            % hemisphere
        end
     
    end    
    
    ieeg_viewLight(v_d(1),v_d(2))
     
    if isequal(mode{m},'SPESclin')
        figureName = fullfile(myDataPath.CCEPpath,'MRI-render','N1_latency_left_CLIN'); 
        title('N1-latency left CLIN')
    
    elseif isequal(mode{m},'SPESprop')
        figureName = fullfile(myDataPath.CCEPpath,'MRI-render','N1_latency_left_PROP'); 
        title('N1-latency left prop')
        
    end
    
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',figureName)
    
    hold off
end


%% RIGHT
% make table with electrodes and the average latency for that electrode
v_d = [96 6];

% make sure electrodes pop out
a_offset = .1*max(abs(allmni_coords(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els = allmni_coords+repmat(a_offset,size(allmni_coords,1),1);      
new_els = els(~isnan(els(:,1)),:);

for m = 1:size(mode,2)

    figure
    gr.faces = Rmnipial_face+1;
    gr.vertices = Rmnipial_vert;
    gr = gifti(gr);
    tH = ieeg_RenderGifti(gr); %#ok<NASGU>

    for pat = 1:size(dataBase,2)
        
        clin_colm = 2*pat-1;                      
        prop_colm = 2*pat; 
    
        if isequal(mode{m},'SPESclin')
            ERs_elec = av_lat_elec(:,clin_colm);
        elseif isequal(mode{m},'SPESprop')
            ERs_elec = av_lat_elec(:,prop_colm);
        end
    
        if all(ismember(all_hemi(start_row(pat):start_row(pat+1),:),'R'))      
            % make table with number of N1's and coordinates of all electrodes
            pat_elec_in_els = start_row(pat,:) : start_row(pat,:)+size(dataBase(pat).agreement_parameter.ERs_elecClin,2)-1 ;
            els_with_N1 = [new_els(pat_elec_in_els,:) ,ERs_elec(1:size(dataBase(pat).agreement_parameter.ERs_elecClin,2),:)];
            [~,idx] = sort(els_with_N1(:,4),'descend');       % Rank/sort based on the 4th column, high to low
            els_ranked = els_with_N1(idx,:);                               % Rank coordinates as well based on number of ERs      
            
                               
            % The darkest color indicates the highest latency for each patient.      
            cm = colormap(flipud(hot(size(unique(els_ranked(:,4)),1)+1)));
            cb = colorbar;
            tickLabels = round(min(unique(els_ranked(:,4))): max(unique(els_ranked(:,4)))/(size(cb.Ticks,2)-1) :max(unique(els_ranked(:,4))));
            
            set(cb, 'TickLabels', cellstr(num2str((tickLabels'))));

            unique_color = 1;
            for elec = 1:size(els_ranked,1)
            
                if elec == 1                        % Only the first elec has to start with a unique color
                    ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                    if els_ranked(elec,4) == els_ranked(elec+1,4)
                        % do nothing, unique_color should remain the same
                    elseif els_ranked(elec,4) ~= els_ranked(elec+1,4)
                        unique_color = unique_color + 1; 
                    end
            
                else
                    if els_ranked(elec,4) == els_ranked(elec-1,4) % if the next stimpair has the same number of N1's as the previous, then give the same color
                        ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                    
                    else
                        ieeg_elAdd(els_ranked(elec,1:3), cm(unique_color,:),12)
                        unique_color = unique_color + 1;        % if the the number of N1's is different from the previous elec, then go to next color
                    end
                end
    
                hold on
            
            end
        
        else 
            % Do nothing because electrodes of this patient are on the right
            % hemisphere
        end
     
    end    
    
    ieeg_viewLight(v_d(1),v_d(2))
     
    if isequal(mode{m},'SPESclin')
        figureName = fullfile(myDataPath.CCEPpath,'MRI-render','N1_latency_right_CLIN'); 
        title('N1-latency right CLIN')
    
    elseif isequal(mode{m},'SPESprop')
        figureName = fullfile(myDataPath.CCEPpath,'MRI-render','N1_latency_right_PROP'); 
        title('N1-latency right prop')
        
    end
    
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',figureName)
    
    hold off
end


end