function fig_average_resp(dataBase, myDataPath)

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
    names_remove(i,:) = dataBase(loc_remove(i)).ccep_clin.sub_label;
end

% Find loc_remove in files 
names = {files.name};

for i = 1:size(names_remove,1)
    files_remove(i,:) = find(contains(names, names_remove(i,:)));
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
        dataBase_clin(R) = load(filename);
        

    elseif isequal(names{R}(strfind(names{R},'task-')+5:strfind(names{R},'task-')+12), 'SPESprop')
        filename = fullfile(myDataPath.dataPath,'derivatives', '/', 'CCEPs/' ,[names{R}]);
        dataBase_prop(R) = load(filename);


    end

end

% /Fridge/users/sifra/derivatives/PRIOSstudy/CCEP_files_allPat
% dit ook nog inladen want hierin zit de ccep_checked met daarin dus of er
% wel of geen ER is.

mode = {'clin','prop'};

data_all_clin = struct;
data_all_prop = struct;


% to indicate the columns in the struct in which the data of all
% patients is saved
pat_c = 1;
pat_p = 1;

for pat = 1:2:size(names,2)
  
    subj_name = names{pat}(strfind(names{pat},'sub-'):strfind(names{pat},'sub-')+10);
    % Load CCEP files to be able to only average the signals that showed an ER
    FilesInFolder = dir(fullfile(myDataPath.CCEP_allpat,[subj_name,'_*']));

%     for m = 1:size(mode,2)        

%         if isequal(mode{m},'clin')
            idx_file_c = contains({FilesInFolder.name},'SPESclin')';

%         elseif isequal(mode{m},'prop')
            idx_file_p = contains({FilesInFolder.name},'SPESprop')';

%         end

        % Load CCEP file
        ccep_clin = load([FilesInFolder(idx_file_c).folder,'/',FilesInFolder(idx_file_c).name],'ccep') ;           
        ccep_prop = load([FilesInFolder(idx_file_p).folder,'/',FilesInFolder(idx_file_p).name],'ccep') ;

        ccep_clin = ccep_clin.ccep;
        ccep_prop = ccep_prop.ccep;

        clin = ccep_clin.n1_peak_sample_check ;
        prop = ccep_prop.n1_peak_sample_check  ;
        fs = 2048;
        ts = 1/fs; 
        clin = ((clin*ts)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
        prop = ((prop*ts)-2)*1000;   

        
        % Only use the ERs that were visually checked.
%         ER_det = find(~isnan(ccep.n1_peak_amplitude_check));

        % For each ER, save the original cc_epoch_sorted REREFERENCED data. 
%         for i = 1:size(ER_det,1)
%             stim = ceil(ER_det(i)/size(ccep.n1_peak_amplitude_check,1));
%             elec = ER_det(i) - size(ccep.n1_peak_amplitude_check,1) * (stim-1);
    
            %%% VOLGENS MIJ KAN IK HIER TUSSENVOEGEN DAT CLIN EN PROP
            %%% ALLEBEI EEN N1-PEAK MOETEN HEBBEN EN DAN ALLEEN DIE DATA
            %%% OPSLAAN

            %%% HOPELIJK KLOPT DE MEDIAN DAN OOK BETER MET DE N1-PEAK
            %%% LATENTIE VAN HET FIGUUR.

            %%% FIGUUR DUS ONDER ELKAAR PLOTTEN EN NIET DOOR ELKAAR HEEN EN
            %%% DAN MET EEN STIPPELLIJNTJE OFZO TER HOOGTE VAN DE N1-PEAK
            %%% VAN DE ' ANDERE' 
        
            i = 1;

      for stimp = 1:size(ccep_prop.stimsets_avg,1)                          % For each stimpair
        for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode

            
        % When both clinical SPES AND propofol SPES show an ER
          if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 

              % Save response data of each response that occurs in clin and prop
              % per patient. DEPTH ELECTRODES in response and stimpair are excluded.
              elec1_stimp = ccep_prop.stimsets_avg(stimp,1);
              elec2_stimp = ccep_prop.stimsets_avg(stimp,2);

              if ~isequal(dataBase_clin(pat).dataBase_clin.tb_channels.group{elec,:}, 'depth') && ~isequal(dataBase_clin(pat).dataBase_clin.tb_channels.group{elec1_stimp,:}, 'depth') && ~isequal(dataBase_clin(pat).dataBase_clin.tb_channels.group{elec2_stimp,:}, 'depth') 
                             
                  data_all_clin(pat).data(i,:) = squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stimp,:))';
                  data_all_prop(pat).data(i,:) = squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stimp,:))';
                  i = i+1;  

              else
                  % Do nothing because response electrode or (one of the)
                  % stimulation pair electrodes is a depth electrode

              end

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

tt = ccep_clin.tt;


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
plot(tt(4147), -291.633,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',6)
% plot(tt(4160), -280.262,'o','MarkerFaceColor',[252/255, 96/255, 57/255],'MarkerEdgeColor',[252/255, 96/255, 57/255],'MarkerSize',6)


% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_clin, fliplr(max_clin)],[17/255, 145/255, 250/255],'EdgeAlpha',0,'FaceAlpha',0.2);

xlim([-0.0140 0.1404])
ylim([-610 216])
xlabel('Time(s)')
ylabel('Potential (\muV)')
% title('Averaged response of all patients during clinical-SPES')
legend('Artefact period','','','Clinical-SPES','N1-peak Clinical-SPES','Location','southeast')


% Create textarrow
annotation(gcf,'textarrow',[0.465503246753247 0.380692640692641],...
    [0.713095238095241 0.713095238095241],'TextBackgroundColor',[1 1 1],...
    'String',{'22.0 ms'},...
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
plot(tt(4155), -288.262,'o','MarkerFaceColor',[252/255, 96/255, 57/255],'MarkerEdgeColor',[252/255, 96/255, 57/255],'MarkerSize',6)
% plot(tt(4147), -252.633,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',6)

% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_prop, fliplr(max_prop)],[252/255, 96/255, 57/255],'EdgeAlpha',0,'FaceAlpha',0.2);

xlim([-0.0140 0.1404])
ylim([-610 216])
xlabel('Time(s)')
ylabel('Potential (\muV)')
sgtitle('Averaged responses of all patients during Clinical-SPES and during Propofol-SPES')
legend('Artefact period','','','Propofol-SPES','N1-peak Propofol-SPES','Location','southeast')


% Create textarrow
annotation(gcf,'textarrow',[0.478354978354978 0.393544372294373],...
    [0.238095238095241 0.238095238095241],'TextBackgroundColor',[1 1 1],...
    'String',{'26.4 ms'},...
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




%% Display clinical response and propofol response in subplots.
% In this example PRIOS03 is used.
pat = 5;
stimp = 1;
elec = 17;

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

plot(tt, squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref(elec,stimp,:,:)),':','Color',[252/255, 96/255, 57/255])
plot(tt, squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'Color',[17/255, 145/255, 250/255],'linewidth',2);   
plot(tt(4153), -563  ,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',4)

xlim([-0.0140 0.1404])
ylim([-1000 500])
xlabel('Time (s)')
ylabel('Potential (\muV)')
%title(sprintf('%s, stimulation pair %s, response electrode %s, Clinical-SPES', subj_name, dataBase_clin(5).dataBase_clin.stimpnames_avg{stimp}, dataBase_clin(5).dataBase_clin.ch{chan}))
title('Clinical-SPES')
legend('Artefact period','','','Separate responses','','','','','','','','','','Averaged signal','N1-peak')

% Create textarrow
annotation(gcf,'textarrow',[0.401222943722944 0.371222943722944],...
    [0.684523809523811 0.684523809523811],'TextBackgroundColor',[1 1 1],...
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

plot(tt, squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref(elec,stimp,[2 7],:)),':','Color',[252/255, 96/255, 57/255])
plot(tt, squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stimp,:)),'Color',[17/255, 145/255, 250/255],'linewidth',2);   
plot(tt(4165), -311,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',4)

xlim([-0.0140 0.1404])
ylim([-700 800])
xlabel('Time (s)')
ylabel('Potential (\muV)')
title('Propofol-SPES')
legend('Artefact period','','','Separate responses','','Averaged signal','N1-peak')

% Create textarrow
annotation(gcf,'textarrow',[0.43858225108225 0.415449134199134],...
    [0.200952380952382 0.200952380952383],'TextBackgroundColor',[1 1 1],...
    'String',{'33.7 ms'},...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');


% Save figure
outlabel=sprintf('clin andprop response PRIOS03_vertical.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

%%
%%% NU PLOT IK ECHT ALLE RESPONSES!!!
%%% IK MOET ALLEEN DE RESPONSES PLOTTEN DIE EEN ER LATEN ZIEN
% clin
data_clin = [mean(dataBase_clin(1).dataBase_clin.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_clin(3).dataBase_clin.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_clin(5).dataBase_clin.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_clin(7).dataBase_clin.cc_epoch_sorted_select_reref_avg),...
    mean(dataBase_clin(9).dataBase_clin.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_clin(11).dataBase_clin.cc_epoch_sorted_select_reref_avg),...
    mean(dataBase_clin(13).dataBase_clin.cc_epoch_sorted_select_reref_avg)];

mean_clin = squeeze(mean(data_clin));
minmax_clin = squeeze(prctile(data_clin),[25,75]);

% prop
data_prop = [mean(dataBase_prop(4).dataBase_prop.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_prop(6).dataBase_prop.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_prop(8).dataBase_prop.cc_epoch_sorted_select_reref_avg),...
    mean(dataBase_prop(10).dataBase_prop.cc_epoch_sorted_select_reref_avg), ...
    mean(dataBase_prop(12).dataBase_prop.cc_epoch_sorted_select_reref_avg),...
    mean(dataBase_prop(14).dataBase_prop.cc_epoch_sorted_select_reref_avg)];

%%%% PRIOS01 HEEFT ER NOG HEEL VEEL RUIS IN ZITTEN TIJDENS PROPOFOL-SPES.
%%%% EVEN CHECKEN OF WEL DE GOEIE REREF OPGESLAGEN IS.
mean_prop = squeeze(mean(data_prop));
minmax_prop = squeeze(prctile(data_prop),[25,50, 75]);

tt = dataBase_clin(1).dataBase_clin.tt;

H=figure();
H.Units = 'normalized';
H.Position = [0.14,0.0625,0.77,0.7];

plot(tt,mean_clin,'r','linewidth',2);  % plot the rereference signal in a solid line
plot(tt,minmax_clin(1),[1,0,0,0.6],'linewidth',2);  % plot the rereference signal in a solid line
plot(tt,minmax_clin(2),[1,0,0,0.6],'linewidth',2);  % plot the rereference signal in a solid line

hold on
plot(tt,mean_prop,'b','linewidth',2);  % plot the rereference signal in a solid line
plot(tt,minmax_prop(1),[0,0,1,0.6],'linewidth',2);  % plot the rereference signal in a solid line
plot(tt,minmax_prop(2),[0,0,1,0.6],'linewidth',2);  % plot the rereference signal in a solid line


xlim([-0.02 0.15])
ylim([-350 310])
xlabel('Time(s)')
ylabel('Potential (\muV)')
title('Averaged response of all patients during clinical-SPES and during propofol-SPES')
legend('Clinical-SPES','Propofol-SPES')

patch([0 0.009 0.009 0],[-800 -800 750 750],[0.6,0.2,0.2],'EdgeAlpha',0)
alpha(0.2)

%% Burst suppression example
% PRios06
% chan is 32
% stimp = 20

H=figure();
H.Units = 'normalized';
H.Position = [0.14,0.0625,0.77,0.7];

subplot(2,1,1)
% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-4000 -4000 7500 7500],[0.9,0.9,0.9],'EdgeAlpha',0)
hold on

plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,[1,5],:)),'Color',[17/255, 145/255, 250/255],'linewidth',2);                % plot the 10 separate stimulations
plot(tt(4152), -1163  ,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',4)
plot(tt(4155), -1070  ,'o','MarkerFaceColor',[17/255, 145/255, 250/255],'MarkerEdgeColor',[17/255, 145/255, 250/255],'MarkerSize',4)

xlim([-0.05 0.15])
ylim([-2000 2000])
xlabel('Time(s)')
ylabel('Potential (\muV)')
legend('Artefact period','Before burst suppression','','N1-peaks','Location','southeast')
title('Before burst suppression')


% Create textarrow
annotation(gcf,'textarrow',[0.487126623376623 0.457126623376623],...
    [0.65357142857143 0.65357142857143],'TextBackgroundColor',[1 1 1],...
    'String',{'N1-peaks before burst suppression'},...
    'LineWidth',1,...
    'FontSize',11,...
    'FontName','FreeSans');



hold off

subplot(2,1,2)
% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-4000 -4000 7500 7500],[0.9,0.9,0.9],'EdgeAlpha',0)
hold on

plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,[2,6],:)),'Color',[252/255, 96/255, 57/255],'linewidth',2);                % plot the 10 separate stimulations
% title(sprintf('Electrode %s, stimulating %s',dataBase.ch{chan},dataBase.stimpnames_avg{stimp}))
title('During burst suppression')
xlim([-0.05 0.15])
ylim([-2000 2000])
legend('Artefact period','During burst suppression','Location','southeast')
xlabel('Time(s)')
ylabel('Potential (\muV)')

sgtitle('Difference in response caused by burst suppression')

% Save figure
outlabel=sprintf('BurstSup_Example.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

end