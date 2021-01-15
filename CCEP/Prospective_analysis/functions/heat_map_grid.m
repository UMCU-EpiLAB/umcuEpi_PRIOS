function heat_map_grid(myDataPath, ccep_clin, agreement_parameter)
% Use heatmaps to display the values of the network characteristics in the
% grid-structure of a patient.
% This can later be used to layer over a MRI to determine whether
% electrodes are on the same gyrus.

% Heatmaps can be smooth or with hard lines, uncomment the wanted method.

subj = extractBetween(ccep_clin.dataName,'sub-','/ses');

if exist(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xls']),'Sheet','matlabsjabloon');
end

% localize electrodes in grid
x = NaN(size(ccep_clin.ch)); 
y = NaN(size(ccep_clin.ch));
elecmat = NaN(size(elec));
% topo=struct;

for i=1:size(elec,1)
    for j=1:size(elec,2)
        if ~ismissing(elec{i,j})
            letter = regexp(elec{i,j},'[a-z,A-Z]');
            number = regexp(elec{i,j},'[1-9]');
            test1 = elec{i,j}([letter,number:end]);
            test2 = [elec{i,j}(letter),'0',elec{i,j}(number:end)];
            if sum(strcmp(ccep_clin.ch,test1))==1
                elecmat(i,j) = find(strcmp(ccep_clin(1).ch,test1));
                y(strcmp(ccep_clin.ch,test1),1) = i;
                x(strcmp(ccep_clin.ch,test1),1)= j;
            elseif sum(strcmp(ccep_clin.ch,test2))==1
                elecmat(i,j) = find(strcmp(ccep_clin.ch,test2));
                y(strcmp(ccep_clin.ch,test2),1) = i;
                x(strcmp(ccep_clin.ch,test2),1)= j;
            else
                error('Electrode is not found')
            end
        end
    end
end

% topo.x = x;
% topo.y = y;
%% For every network characteristic
% For prop and clin
mode = {'Indegree','Outdegree','BC'};

for J = 1:size(mode,2) 
    
    figure('Position',[284,4,1309,1052]);
    
    if strcmp(mode{J},'Indegree')
        parclin(:,1) = ccep_clin.ch;
        parclin(:,2) = num2cell(agreement_parameter.indegreeN_Clin)';
        
        parprop(:,1) = ccep_clin.ch;
        parprop(:,2) = num2cell(agreement_parameter.indegreeN_Prop)';
                
    elseif strcmp(mode{J},'BC')
        parclin(:,1) = ccep_clin.ch;
        parclin(:,2) = num2cell(agreement_parameter.BCN_Clin)';
        
        parprop(:,1) = ccep_clin.ch;
        parprop(:,2) = num2cell(agreement_parameter.BCN_Prop)';
                
    elseif strcmp(mode{J},'Outdegree')
        parclin(:,1) = ccep_clin.ch;
        parclin(:,2) = num2cell(agreement_parameter.outdegreeN_Clin)';
        
        parprop(:,1) = ccep_clin.ch;
        parprop(:,2) = num2cell(agreement_parameter.outdegreeN_Prop)';
        
    end
      
%%%%%%%%%%%%%%%%% Clinical SPES %%%%%%%%%%%%%%%%%%%%%%%%
val_mat = elec;

% Match the electrode names with the value of the network characteristic
for rows = 1:size(elec,1)
    
    for col = 1:size(elec,2)
    
        if ~isnan(elecmat(rows,col))
           elec_loc = find(ismember(parclin(:,1), elec(rows,col)));                % Determine the indegree value of the elec
           val = parclin{elec_loc,2}; %#ok<FNDSB>
           val_mat(rows,col) = {val};          
    
        end  
    end 
end

% Replace missing with zero to later transform the cell matrix to double
% matrix
mask = cellfun(@ismissing, val_mat);
val_mat(mask) = {0};    

% If sum of row is zero, remove row, same for column (to reduce space)
for rows = size(val_mat,1):-1:1
    if max([val_mat{rows,:}]) == 0
        val_mat(rows,:) = [];
    end
end

for cols = size(val_mat,2):-1:1
    if max([val_mat{:,cols}]) == 0
        val_mat(:,cols) = [];
    end
end

% Transform to a double matrix to perform imagesc with colormap
mat = flip(cell2mat(val_mat));
 
% Axes of clinical-SPES
subplot(2,1,1)
pcolor(mat)         % Smooth transition between electrodes
shading interp
%imagesc(mat)        % Clear transision between electrodes
colormap(flipud(hot))               % Darker colors is higher, light is low.
colorbar
axes1 = gca;
axes1.YTick = [];                                               % Remove numbers on y-axis
axes1.XTick = [];                                               % Remove numbers on x-axis
axes1.XColor = 'none';                                          % Remove line indicating the x-axis
axes1.YColor = 'none';                                          % Remove line indicating the y-axis
title(sprintf('Clinical-SPES, %s, %s',mode{J},subj{:}))

%%%%%%%%%%%%%%%%% Propofol SPES %%%%%%%%%%%%%%%%%%%%%%%%
val_mat_prop = elec;

% Match the electrode names with the value of the network characteristic
for rows = 1:size(elec,1)
    
    for col = 1:size(elec,2)
    
        if ~isnan(elecmat(rows,col))
           elec_loc = find(ismember(parprop(:,1), elec(rows,col)));                % Determine the indegree value of the elec
           val = parprop{elec_loc,2}; %#ok<FNDSB>
           val_mat_prop(rows,col) = {val};          
    
        end  
    end 
end

% Replace missing with zero to later transform the cell matrix to double
% matrix
mask = cellfun(@ismissing, val_mat_prop);
val_mat_prop(mask) = {0};    

% If sum of row is zero, remove row, same for column (to reduce space)
for rows = size(val_mat_prop,1):-1:1
    if max([val_mat_prop{rows,:}]) == 0
        val_mat_prop(rows,:) = [];
    end
end

for cols = size(val_mat_prop,2):-1:1
    if max([val_mat_prop{:,cols}]) == 0
        val_mat_prop(:,cols) = [];
    end
end

% Transform to a double matrix to perform imagesc with colormap
% Flip to make sure the orientation is as expected
mat = flip(cell2mat(val_mat_prop));
 
% Axes of Propofol-SPES
subplot(2,1,2)
pcolor(mat)         % Smooth transition between electrodes
shading interp
% imagesc(mat)        % Clear transision between electrodes
colormap(flipud(hot))               % Darker colors is higher, light is low.
colorbar
axes2 = gca;
axes2.YTick = [];                                               % Remove numbers on y-axis
axes2.XTick = [];                                               % Remove numbers on x-axis
axes2.XColor = 'none';                                          % Remove line indicating the x-axis
axes2.YColor = 'none';                                          % Remove line indicating the y-axis
title(sprintf('Propofol-SPES, %s, %s',mode{J},subj{:}))

% Save figure
outlabel=sprintf('sub-%s_%s.jpg',...
    subj{1},mode{J});
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/HeatMap_grid/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')




%% MRI figures with a layer of the network characteristics
% Currently only for PRIOS06 becuase this patient had the 'easiest' MRI.

if ismember('PRIOS06', subj)
    
    figure()
    MRI = imread('PRIOS06.jpg');
    imshow(MRI)
    title(sprintf('%s Indegree SPES-clin',subj{:}))

    % Get the points of the elektrodes in the CT scan
    % Make sure to select the electrodes in the same order as the electrode
    % names in dataBase.ch!!! (check excel sjabloon when nessecary)
%     [xi,yi] = getpts;                                                      % When done, press enter twice to save the coordinates
    hold on
    colormap(flipud(hot))  
    im = scatter(xi, yi, 260, agreement_parameter.indegreeN_Clin(:,1:size(xi,1))' ,'filled');
    alpha(im,0.5);

    % Save figure
    outlabel=sprintf('sub-%s_SPES_Clin Indegree.jpg',...
        subj{1});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/HeatMap_grid/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')


end
end