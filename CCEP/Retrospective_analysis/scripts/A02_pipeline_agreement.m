%% Script pipeline_preproces.m should be performed first to obtain the correct documents.
% This script analyses the agreement when 2 or 10 stimuli are averaged and
% N1s are detected.

%%
clear; 
% Select all patients
cfg.sub_labels = {'sub-RESP0701','sub-RESP0702','sub-RESP0703','sub-RESP0706','sub-RESP0724','sub-RESP0728'}; 

% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% Load all ccep files in the folder CCEP_files_allPat
% Be aware that for some patients, SPES is saved in multiple runs, but
% these are already merged in the previous step

files = dir(fullfile(myDataPath.CCEP_allpat));

% Create database with the CCEP information of all patients of all runs and
% both protocols (2 and 10 stims)
dataBase = struct;                      
for i=1:size(cfg.sub_labels,2)                                                % DataBase must have the size of the number of runs 
    dataBase(i).sub_label = cfg.sub_labels{i};        
    
    % Find rows with the sub_label of interest 
    respLoc = find(contains({files(:).name},cfg.sub_labels{i}));              
    
    % load all both the 10 stimuli and 2 stimuli of the patient
    for j=1:size(respLoc,2)                                                   % number of rows with the run_label of interest
       if contains(files(respLoc(j)).name,'10stims') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep10 = ccep10;
          dataBase(i).filename10 = files(respLoc(j)).name;
          
       elseif contains(files(respLoc(j)).name,'2stims') 
          load(fullfile(files(respLoc(j)).folder,files(respLoc(j)).name));
          dataBase(i).ccep2 = ccep2;   
          dataBase(i).filename2 = files(respLoc(j)).name;
       end
    end
end

% small cleanup
clear respLoc k j files ccep10 ccep2 

%% Determine the agreement between 2 and 10 stims per run
close 
clc

for subj = 1:size(dataBase,2)
    
    % Find the 2 runs matching.
    runs(1).ccep = dataBase(subj).ccep10;
    runs(1).name = dataBase(subj).filename10;                                                  
    runs(1).sub_label = dataBase(subj).sub_label;
    runs(2).ccep = dataBase(subj).ccep2;
    runs(2).name = dataBase(subj).filename2;
    runs(2).sub_label = dataBase(subj).sub_label;
    
    % Determine the agreement between the two matching runs
    agreement = determine_agreement(runs);          
    
    dataBase(subj).agreement = agreement;
    
    fprintf('%s Overall agreement = %1.2f, positive agreement = %1.2f, negative agreement = %1.2f \n',...
        dataBase(subj).sub_label, agreement.agreement_run.OA, agreement.agreement_run.PA, agreement.agreement_run.NA)
end

% short clean up
clear subj agreement runs


%% Determine the location of the ones (ER vs. No-ER)
% ccep10 is only used for the channels and stimpairs and those are
% equal for 2 and 10 stimuli so does not matter which database is used.

for subj = 1:size(dataBase,2)
    ccep10 = dataBase(subj).ccep10;
    agreement = dataBase(subj).agreement;
    
    LocOnes = find_ones(ccep10, agreement.agreement_run);
    dataBase(subj).LocOnes = LocOnes;
end

% clean up
clear ccep10 agreement LocOnes subj


%% Rewrite the adjacency matrices to electrode x electrode instead of electrode x stimulation pair
% This is necessary to calculate the indegree, outdegree and BC.
close all;

for subj = 1:size(dataBase,2)
    
    dataBase(subj).ccep2 = rewrite_Amat(dataBase(subj).ccep2, dataBase(subj).agreement.Amat2);
    dataBase(subj).ccep10 = rewrite_Amat(dataBase(subj).ccep10, dataBase(subj).agreement.Amat10);
end

%% Determine the indegree, outdegree, Betweenness centrality, the number of ERs per stimpair and the number of ERs per electrode
% The variables are saved in an excel in the run folder of the subject number
close all;
clc

for subj = 1:size(dataBase,2)
    dataBase(subj).agreement_parameter = agreement_parameters(dataBase(subj).agreement, ...
        dataBase(subj).ccep2, dataBase(subj).ccep10);
    
    dataBase(subj).statistics = statistical_agreement(myDataPath, dataBase(subj).agreement_parameter, dataBase(subj).ccep10);
end


%% Scatter plot of the network parameters
% function corrplot cannot be used because I don't have the Econometrics
% Toolbox, though refline also works perfectly.
close all
scatter_networkPar(dataBase, myDataPath)

%% Load electrodes positions (xlsx/electrodes.tsv)
plot_fig = 'n';                 % 'n' when you don't want all ER responses per stim, and the SOZ to be plot
close all;

for subj = 1:size(dataBase,2)
    visualise_gridstructure(myDataPath, dataBase(subj).ccep10, dataBase(subj).ccep2, dataBase(subj).agreement_parameter,plot_fig);
end

%% Display agreement parameters on brain image

%%% DIT IS EVEN VOOR HET VOORBEELD VOOR IN MIJN VERSLAG!!!
%%% waarbij ik pijlen teken tussen electrodes om de indegree, outdegree en
%%% BC uit te leggen. Maar ben niet tevreden over het resultaat.

close all;
clear coordinates

% Show CT scan with the grid (patient specific)
[I] = imread('CT met grid.jpg');            
figure1 = figure();
fig = imshow(I);
set(figure1, 'Position', [319,7,1241,935]);        % enlarge figure
hold on

% Get the points of the elektrodes in the CT scan
% Make sure to select the electrodes in the same order as the electrode
% names in dataBase.ch!!! (check excel sjabloon when nessecary)
[xi,yi] = getpts;                                                      % When done, press enter twice to save the coordinates
coordinates(:,1:2) = [xi(:),yi(:)];

% Put random electrode names next to the selected points in the CT scan
set(fig, 'AlphaData', 0.6);         % Set transparancy
% array of random EL names

for i = 1:size(xi,1)
    name_ch{i,:} = ['LC' num2str(i,'%2d')];
    
end
for i = 16:48
    name_ch{i,:} = ['C' num2str(i-16,'%2d')];
end

text(((xi)),yi, name_ch, 'FontSize',8,'FontWeight','bold')           % Does not matter wheter ccep10 or ccep2 is chosen, ch is the same



% PLot the indegree between the electrodes 
% When the elec_mat shows a 1, there is a connection between the electrodes
% Therefore a line is drawn between the electrodes
temp_mat = dataBase(subj).ccep10.elec_Amat(1:length(xi),1:length(xi));
% Replace all twos in the elec x elec matrix with 1 since this is not of
% interest now
twos = find(ismember(temp_mat,2));
temp_mat(twos) = 1;


temp_mat = zeros(size(dataBase(subj).ccep10.elec_Amat(1:length(xi),1:length(xi))))
% 10 stims
%%% dit nog met 'mode'  proberen in 1 loop te zetten
lineWidth_ind = linspace(min(dataBase(subj).agreement_parameter.indegreeN_10), max(dataBase(subj).agreement_parameter.indegreeN_10), 4);        % devide the number of ERs to fit the lineWidth 
bins_lineWidth_ind= discretize(dataBase(subj).agreement_parameter.indegreeN_10, lineWidth_ind)*2;

% Thickness and transperancy of lines indicates the indegree
for elec1 = 1:length(temp_mat)                                     % This should be the same as the number of channels!!
    for elec2 = 1:length(temp_mat)
       
        if temp_mat(elec1,elec2)== 1                                % When there is a link between elec 1 and elec 2
            
            
            % Lines between electrodes with ER connection 
            x = [xi(elec1),yi(elec1)];
            y = [xi(elec2),yi(elec2)];
            dp = y - x;

            arrow_elec = quiver(x(1),x(2),dp(1),dp(2),0,'Color',[0.05 0.05 0.05],'LineWidth', 1.5,'AutoScaleFactor',06) ;              %(bins_lineWidth_ind(elec1)));
            %arrow_elec.Color(4) = bins_lineWidth_ind(elec1)*0.05;                  %the lower the more transparant
       end
    end
end
