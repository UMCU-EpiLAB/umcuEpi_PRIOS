%% display connections

%% pre-allocation
clear
config_CCEP

%% set paths
cfg.mode = 'retro';
myDataPath = setLocalDataPath(cfg);

%% select run
% choose between available runs
files = dir(fullfile(myDataPath.dataPath,cfg.sub_labels{1}, cfg.ses_label,'ieeg',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_events.tsv']));
names = {files.name};
strings = cellfun(@(x) x(strfind(names{1},'run-'):strfind(names{1},'run-')+9), names, 'UniformOutput', false);
stringsz = [repmat('%s, ',1,size(strings,2)-1),'%s'];

cfg.run_label = {input(sprintf(['Choose one of these runs: \n' stringsz '\n'],strings{:}),'s')}; % Chosen run is in cfg.run_label

clear files names strings stringsz

%% Load all CCEP set 
% Does not matter whether you take the 10 stims or the 2 stims
% Only the channels and stimpairs are used, which should be similar in both
% stimulation protocols.

files = dir(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,cfg.run_label{1}));
n=1; 
runs = cell(1);

for i=1:size(files,1)
    if contains(files(i).name,'CCEP_') && n==1
        loadfile = load(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,cfg.run_label{1},[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',cfg.run_label{:},'_CCEP_2stims.mat']));
        ccep = loadfile.ccep2;
        runs{n} = files(i).name;
        n=n+1;
    elseif contains(files(i).name ,'CCEP_') && n>1
        loadfile = load(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,cfg.run_label{1},[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',cfg.run_label{:},'_CCEP_2stims.mat']));
        ccep(n) = loadfile.ccep2;
        runs{n} = files(i).name;
        n=n+1;
    end
end


%% load electrodes positions (xlsx/electrodes.tsv)

subj = cfg.sub_labels{1}(5:end);

if exist(fullfile(myDataPath.elec_input,[subj,'_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj,'_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj,'_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj,'_elektroden.xls']),'Sheet','matlabsjabloon');
end

% localize electrodes in grid
x = NaN(size(ccep(1).ch)); 
y = NaN(size(ccep(1).ch));
elecmat = NaN(size(elec));
topo=struct;

for i=1:size(elec,1)
    for j=1:size(elec,2)
        if ~ismissing(elec{i,j})
            letter = regexp(elec{i,j},'[a-z,A-Z]');
            number = regexp(elec{i,j},'[1-9]');
            test1 = elec{i,j}([letter,number:end]);
            test2 = [elec{i,j}(letter),'0',elec{i,j}(number:end)];
            if sum(strcmp(ccep(1).ch,test1))==1
                elecmat(i,j) = find(strcmp(ccep(1).ch,test1));
                y(strcmp(ccep(1).ch,test1),1) = i;
                x(strcmp(ccep(1).ch,test1),1)= j;
            elseif sum(strcmp(ccep(1).ch,test2))==1
                elecmat(i,j) = find(strcmp(ccep(1).ch,test2));
                y(strcmp(ccep(1).ch,test2),1) = i;
                x(strcmp(ccep(1).ch,test2),1)= j;
            else
                error('Electrode is not found')
            end
        end
    end
end

topo.x =x;
topo.y=y;


%% Plot electrodes which respond differently to the 2 stimuli and the 10 stimuli protocol

for stimp = 1:size(ccep(1).stimsets_avg)                   % Number of stimulation pairs (columns)
    stimnum = ccep(1).stimsets_avg(stimp,:);            % Stimulation pair numbers for column number (stimp)
    resp = compare_mat(:,stimp);                    % matrix with one and zero for ER and non ER, respectively

        
        figure(1),
        % plot all electrodes
        plot(topo.x,topo.y,'ok','MarkerSize',15)
        hold on
        
        
        % plot electrodes showing CCEPs in green (CCEP = 2 because 2 and 10 stimulations are compared)
        chan = find(resp==2);
        plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
            'MarkerFaceColor','g','MarkerEdgeColor','k')
        
        % plot electrodes showing CCEPs in green (CCEP = 2 because 2 and 10 stimulations are compared)
        chan = find(resp==0);
        plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
            'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','k')
                
        % plot electrodes showing CCEPs in one of the two stimulations in red 
        chan = find(resp==1);
        plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
            'MarkerFaceColor','r','MarkerEdgeColor','k')
       
        % plot stimulation pair in yellow
        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor','y','MarkerEdgeColor','k')
        end
        
        hold off
        
        % add electrode names
        text(topo.x,topo.y,ccep(1).ch)
        
        ax = gca;
        xlim([min(topo.x)-2, max(topo.x)+2])
        ylim([min(topo.y)-2, max(topo.y)+2])
        title(sprintf('CCEP responses after stimulating %s-%s', ccep(1).ch{stimnum(1)}, ccep(1).ch{stimnum(2)}))
        
        ax.YDir = 'reverse';
        ax.YTick = [];
        ax.XTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.Units = 'normalized';
        ax.Position = [0.1 0.1 0.8 0.8];
        outlabel=sprintf('Stimpair%s-%s.jpg',...
            ccep(1).ch{stimnum(1)},ccep(1).ch{stimnum(2)});
        
        path = fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label, cfg.run_label{:});
        if ~exist([path,'/figures/'], 'dir')
            mkdir([path,'/figures/']);
        end
        
        saveas(gcf,[path,'/figures/',outlabel],'jpg')
        
end
