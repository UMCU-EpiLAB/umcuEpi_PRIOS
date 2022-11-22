function fig_average_resp(dataBase, myDataPath)
% Create a figure of the median of all responses of all patients.

files = dir(fullfile(myDataPath.dataPath, 'derivatives', '/', 'CCEPs' ,'/','*_CCEP.mat'));

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

loc_remove = find(dataBase_remove == 1);

for i = 1:size(loc_remove,2)    
    names_remove(i,:) = dataBase(loc_remove(i)).ccep_clin.sub_label; %#ok<AGROW> 
end

% Find loc_remove in files 
names = {files.name};

for i = 1:size(names_remove,1)
    files_remove(i,:) = find(contains(names, names_remove(i,:))); %#ok<AGROW> 
end

remove_files = unique(files_remove);

dataBase(dataBase_remove ==1) = [];
files(remove_files',:) = [];

names = {files.name};
strings = cell(size(names));

% ONly select the subjects of whom a clinical-SPES and propofol-SPES is
% available
for n = 1:size(names,2)
    strings{n} = names{n}(strfind(names{n},'sub-'):strfind(names{n},'sub-')+10);
end

for R = 1:size(strings,2)
    if isequal(names{R}(strfind(names{R},'task-')+5:strfind(names{R},'task-')+12), 'SPESclin')
        filename = fullfile(myDataPath.dataPath,'derivatives', '/', 'CCEPs/' ,[names{R}]);
        dataBase_clin(R) = load(filename); %#ok<AGROW> 
        
    elseif isequal(names{R}(strfind(names{R},'task-')+5:strfind(names{R},'task-')+12), 'SPESprop')
        filename = fullfile(myDataPath.dataPath,'derivatives', '/', 'CCEPs/' ,[names{R}]);
        dataBase_prop(R) = load(filename); %#ok<AGROW> 
    end
end

% Save cc_epoch_sorted_select_reref_avg WITH A CHECKED CCEP of all patients 
% in one file to be averaged and plot in 1 figure. 
data_all_clin = struct;
data_all_prop = struct;

for pat = 1:2:size(names,2)
     
    % Use N1-checked files to find the responses with an N1.
    clin = dataBase(round(pat/2)).ccep_clin.n1_peak_sample;
    prop = dataBase(round(pat/2)).ccep_prop.n1_peak_sample;
    fs = 2048;
    ts = 1/fs; 
    clin = ((clin*ts)-2)*1000;                                            % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
    prop = ((prop*ts)-2)*1000;   
        
    i = 1;

    ccep_clin = dataBase(round(pat/2)).ccep_clin;
    ccep_prop = dataBase(round(pat/2)).ccep_prop;

    for stimp = 1:size(ccep_clin.stimsets_avg,1)                          % For each stimpair
        for elec = 1:size(ccep_prop.ch,1)                                 % For each electrode
            
        % Make sure that the electrodes of the stimulation pair are
        % neither bad or depth electrodes
        elec1_stimp = dataBase(round(pat/2)).ccep_prop.stimsets_avg(stimp,1);
        elec2_stimp = dataBase(round(pat/2)).ccep_prop.stimsets_avg(stimp,2);

        % When both clinical SPES AND propofol SPES show an ER
          if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) && isequal(ccep_clin.tb_channels.status{elec}, 'good') && ~isequal(ccep_clin.tb_channels.group{elec1_stimp,:}, 'depth') && ~isequal(ccep_clin.tb_channels.group{elec2_stimp,:}, 'depth') 

             % Responses on DEPTH and BAD electrodes are already changed to NAN
             data_all_clin(pat).data(i,:) = squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stimp,:))';
             data_all_prop(pat).data(i,:) = squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stimp,:))';
             i = i+1;  

          end 
        end
    end
end

data_clin = prctile(vertcat(data_all_clin(:).data),50);
min_clin= prctile(vertcat(data_all_clin(:).data),25);
max_clin = prctile(vertcat(data_all_clin(:).data),75);


data_prop = prctile(vertcat(data_all_prop(:).data),50);
min_prop = prctile(vertcat(data_all_prop(:).data),25);
max_prop = prctile(vertcat(data_all_prop(:).data),75);

tt = dataBase_clin(1).dataBase_clin.tt;

% Create figure with the average of all signals with an N1-peak
H=figure();
H.Units = 'normalized';
H.Position = [0.14,0.0625,0.77,0.7];

% Create a patch in which the stimulation artefact is shown
subplot(2,1,1)

patch([0 0.009 0.009 0],[-800 -800 750 750],[0.9,0.9,0.9],'EdgeAlpha',0,'FaceAlpha',0.5)
hold on
vline(0.009,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);
vline(0,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);

plot(tt,data_clin,'color',[17/255, 145/255, 250/255],'linewidth',2);  % plot the rereference signal in a solid line
hold on
plot(tt(4149), -314.324,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',6)

% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_clin, fliplr(max_clin)],[17/255, 145/255, 250/255],'EdgeAlpha',0,'FaceAlpha',0.2);

xlim([-0.0140 0.1404])
ylim([-610 216])
xlabel('Time(s)')
ylabel('Potential (\muV)')
legend('Artefact period','','','Clinical-SPES','N1-peak Clinical-SPES','Location','southeast')

% Create a new pair of axes inside current figure
axes('position',[.55 .620 .17 .11])
box on % put box around new pair of axes
indexOfInterest = (tt < tt(4159)) & (tt > tt(4139)); % range of t near peak
plot(tt(indexOfInterest),data_clin(indexOfInterest),'color',[17/255, 145/255, 250/255],'linewidth',2) % plot on new axes
axis tight
ylim([-316.9765625,-301.408735795455])
hold on
plot(tt(4149), -314.324,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',6)

% Create textarrow towards N1-peak
annotation(gcf,'textarrow',[0.525503246753247 0.380692640692641],...
    [0.713095238095241 0.713095238095241],'TextBackgroundColor',[1 1 1],...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');

% Prop
subplot(2,1,2)
% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-800 -800 750 750],[0.9,0.9,0.9],'EdgeAlpha',0,'FaceAlpha',0.5)
hold on
vline(0.009,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);
vline(0,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);

plot(tt,data_prop,'color',[252/255, 96/255, 57/255],'linewidth',2);  % plot the rereference signal in a solid line
hold on
plot(tt(4161), -308.262,'o','MarkerFaceColor',[252/255, 96/255, 57/255],'MarkerEdgeColor',[252/255, 96/255, 57/255],'MarkerSize',6)

% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_prop, fliplr(max_prop)],[252/255, 96/255, 57/255],'EdgeAlpha',0,'FaceAlpha',0.2);

xlim([-0.0140 0.1404])
ylim([-610 216])
xlabel('Time(s)')
ylabel('Potential (\muV)')
sgtitle('Averaged responses of all patients during Clinical-SPES and during Propofol-SPES')
legend('Artefact period','','','Propofol-SPES','N1-peak Propofol-SPES','Location','southeast')

% Create a new pair of axes inside current figure
axes('position',[.55 .15 .17 .11])
box on % put box around new pair of axes
indexOfInterest = (tt < tt(4170)) & (tt > tt(4150)); % range of t near peak
plot(tt(indexOfInterest),data_prop(indexOfInterest),'color',[252/255, 96/255, 57/255],'linewidth',2) % plot on new axes
axis tight
ylim([-312.9765625,-285.408735795455])
hold on
plot(tt(4161), -308.262,'o','MarkerFaceColor',[252/255, 96/255, 57/255],'MarkerEdgeColor',[252/255, 96/255, 57/255],'MarkerSize',6)

% Create textarrow towards N1-peak
annotation(gcf,'textarrow',[0.528354978354978 0.393544372294373],...
    [0.238095238095241 0.238095238095241],'TextBackgroundColor',[1 1 1],...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');

% Save figure
outlabel=sprintf('averagedRESPallPat.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

% savefig(gcf,[path,outlabel])


%% Display ONE example of a clinical and propofol response in subplots.
pat = 1;
stimp = 24;
elec = 36;

tt = dataBase_clin(pat).dataBase_clin.tt;

close all
H=figure(1);
H.Units = 'normalized';
H.Position = [0.227083333333333,0.000833333333333,0.392708333333333,0.706666666666667]; %[0.14,0.0625,0.77,0.7];

subplot(2,1,1)
% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-1000 -1000 750 750],[0.9,0.9,0.9],'EdgeAlpha',0)
hold on
vline(0.009,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);
vline(0,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);

plot(tt, squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref(elec,stimp,:,:)),':','Color',[.5 .5 .5],'LineWidth',1)
plot(tt, squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'Color',[17/255, 145/255, 250/255],'linewidth',2);   
plot(tt(4153), -563  ,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',4)

xlim([-0.0140 0.1404])
ylim([-1000 400])
xlabel('Time (s)')
ylabel('Potential (\muV)')
%title(sprintf('%s, stimulation pair %s, response electrode %s, Clinical-SPES', subj_name, dataBase_clin(5).dataBase_clin.stimpnames_avg{stimp}, dataBase_clin(5).dataBase_clin.ch{chan}))
title('Clinical-SPES')
legend('Artefact period','','','Separate responses','','','','','','','','','','Averaged signal','N1-peak')

% Create textarrow towards N1-peak
annotation(gcf,'textarrow',[0.402549203669894 0.372549203669894],...
    [0.638533243486074 0.638533243486074],'TextBackgroundColor',[1 1 1],...
    'String',{'27.8 ms'},...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');

subplot(2,1,2)
% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-800 -800 1000 1000],[0.9,0.9,0.9],'EdgeAlpha',0)
hold on
vline(0.009,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);
vline(0,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);

plot(tt, squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref(elec,stimp,[2 7],:)),':','Color', [.5 .5 .5],'LineWidth',1.5)
plot(tt, squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'Color',[252/255, 96/255, 57/255],'linewidth',2);   
plot(tt(4165), -311,'o','MarkerFaceColor',[252/255, 96/255, 57/255],'MarkerEdgeColor',[252/255, 96/255, 57/255],'MarkerSize',4)

xlim([-0.0140 0.1404])
ylim([-767 532])
xlabel('Time (s)')
ylabel('Potential (\muV)')
title('Propofol-SPES')
legend('Artefact period','','','Separate responses','','Averaged signal','N1-peak')

% Create textarrow towards N1-peak
annotation(gcf,'textarrow',[0.434603471241401 0.401856763925731],...
    [0.228075022461816 0.227594339622642],'TextBackgroundColor',[1 1 1],...
    'String',{'33.7 ms'},...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');


% Save figure
outlabel=sprintf('clin andprop response_vertical.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

end