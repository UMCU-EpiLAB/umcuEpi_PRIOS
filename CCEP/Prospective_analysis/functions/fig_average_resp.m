function fig_average_resp(myDataPath)

files = dir(fullfile(myDataPath.dataPath, 'derivatives', '/', 'CCEPs' ,'/',['*_CCEP.mat']));

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

    for m = 1:size(mode,2)        

        if isequal(mode{m},'clin')
            idx_file = contains({FilesInFolder.name},'SPESclin')';

        elseif isequal(mode{m},'prop')
            idx_file = contains({FilesInFolder.name},'SPESprop')';

        end

        % Load CCEP file
        load([FilesInFolder(idx_file).folder,'/',FilesInFolder(idx_file).name],'ccep')            
        
        % Only use the ERs that were visually checked.
        ER_det = find(~isnan(ccep.n1_peak_amplitude_check));

        % For each ER, save the original cc_epoch_sorted REREFERENCED data. 
        for i = 1:size(ER_det,1)
            stim = ceil(ER_det(i)/size(ccep.n1_peak_amplitude_check,1));
            elec = ER_det(i) - size(ccep.n1_peak_amplitude_check,1) * (stim-1);
    
            if isequal(mode{m},'clin')
                 data_all_clin(pat_c).data(i,:) = squeeze(dataBase_clin(pat).dataBase_clin.cc_epoch_sorted_select_reref_avg(elec,stim,:))';

            elseif isequal(mode{m},'prop')
                 data_all_prop(pat_p).data(i,:) = squeeze(dataBase_prop(pat+1).dataBase_prop.cc_epoch_sorted_select_reref_avg(elec,stim,:))';
                 
            end

        end 
    
    end
    pat_c = pat_c + 1;
    pat_p = pat_p + 1;
end


data_clin = prctile(vertcat(data_all_clin(:).data),50);
min_clin= prctile(vertcat(data_all_clin(:).data),25);
max_clin = prctile(vertcat(data_all_clin(:).data),75);


data_prop = prctile(vertcat(data_all_prop(:).data),50);
min_prop = prctile(vertcat(data_all_prop(:).data),25);
max_prop = prctile(vertcat(data_all_prop(:).data),75);

tt = ccep.tt;


H=figure();
H.Units = 'normalized';
H.Position = [0.14,0.0625,0.77,0.7];

% Create a patch in which the stimulation artefact is shown
patch([0 0.009 0.009 0],[-800 -800 750 750],[0.9,0.9,0.9],'EdgeAlpha',0,'FaceAlpha',0.5)
hold on
vline(0.009,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);
vline(0,'linewidth',2,'linestyle','--','color',[0.9,0.9,0.9]);


plot(tt,data_clin,'r','linewidth',2);  % plot the rereference signal in a solid line

hold on
plot(tt(4142), -177.7661,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)

% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_clin, fliplr(max_clin)],[0.6,0.2,0.2],'EdgeAlpha',0,'FaceAlpha',0.2);

% Prop
plot(tt,data_prop,'b','linewidth',2);  % plot the rereference signal in a solid line
hold on
plot(tt(4160), -239.6932,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',4)

% Fill between 25% and 75%
patch([tt, fliplr(tt)],[min_prop, fliplr(max_prop)],[0.2,0.2,0.6],'EdgeAlpha',0,'FaceAlpha',0.2);

xlim([-0.0140 0.1404])
ylim([-516 220])
xlabel('Time(s)')
ylabel('Potential (\muV)')
title('Averaged response of all patients during clinical-SPES and during propofol-SPES')
legend('Artefact period','','','Clinical-SPES','','','Propofol-SPES')





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

end