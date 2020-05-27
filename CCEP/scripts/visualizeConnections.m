%% display connections

%% pre-allocation
clear
config_CCEP

%% set paths
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

%% load all CCEP set LOOK AT THE NAME!! 2 VERSUS 10 STIMS!!!

files = dir(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,cfg.run_label{1}));
n=1; 
runs = cell(1);

for i=1:size(files,1)
    if contains(files(i).name ,'run-') && n==1
        loadfile = load(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,cfg.run_label{1},files(i).name,[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',files(i).name,'_CCEP_10stims.mat']));
        ccep = loadfile.ccep;
        runs{n} = files(i).name;
        n=n+1;
    elseif contains(files(i).name ,'run-') && n>1
        loadfile = load(fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,files(i).name,[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',files(i).name,'_CCEP_10stims.mat']));
        ccep(n) = loadfile.ccep;
        runs{n} = files(i).name;
        n=n+1;
    end
end


%% load all CCEP detected by automatic detector

[OA, PA, NA, compare_mat] = determine_agreement(myDataPath,cfg);

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


%% CCEP responses to specific stimulus
% for run = 1:size(ccep,2)
%     adj_matrix{run,1} = ~isnan(ccep(run).n1_peak_sample);
%     ccep(run).adj_matrix = adj_matrix{run,1};
% end

%<<<<<<< HEAD
for stimp = 1:size(compare_mat,2)      % Number of stimulation pairs (columns)
    stimnum = ccep(run).cc_stimsets(stimp,:);       % Stimulation pair numbers for column number (stimp)
    resp = ccep(run).adj_matrix(:,stimp);           % matrix with one and zero for ER and non ER, respectively

    figure(1),
    % plot all electrodes
    plot(topo.x,topo.y,'ok','MarkerSize',15)
    hold on
    % plot stimulation pair in yellow
    for chan=1:2
        plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
            'MarkerFaceColor','y','MarkerEdgeColor','k')

for run = 1:size(ccep,2)
    for stimp = 1:size(ccep(run).checked,2)
        stimnum = ccep(run).cc_stimsets(stimp,:);
        resp = ccep(run).checked(:,stimp);
        
        figure(1),
        % plot all electrodes
        plot(topo.x,topo.y,'ok','MarkerSize',15)
        hold on
        % plot stimulation pair in yellow
        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor','y','MarkerEdgeColor','k')
        end
        
        % plot electrodes showing CCEPs in green
        chan = find(resp==1);
        plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
            'MarkerFaceColor','g','MarkerEdgeColor','k')
        hold off
        
        % add electrode names
        text(topo.x,topo.y,ccep(run).ch)
        
        ax = gca;
        xlim([min(topo.x)-2, max(topo.x)+2])
        ylim([min(topo.y)-2, max(topo.y)+2])
        title(sprintf('CCEP responses after stimulating %s-%s', ccep(run).ch{stimnum(1)}, ccep(run).ch{stimnum(2)}))
        
        ax.YDir = 'reverse';
        ax.YTick = [];
        ax.XTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.Units = 'normalized';
        ax.Position = [0.1 0.1 0.8 0.8];
        outlabel=sprintf('Stimpair%s-%s.jpg',...
            ccep(run).ch{stimnum(1)},ccep(run).ch{stimnum(2)});
        
        path = fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,runs{run});
        if ~exist([path,'/figures/'], 'dir')
            mkdir([path,'/figures/']);
        end
        
        saveas(gcf,[path,'/figures/',outlabel],'jpg')
        
    end

    % plot electrodes showing CCEPs in green
    chan = find(resp==1);
    plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
        'MarkerFaceColor','g','MarkerEdgeColor','k')
    hold off

    % add electrode names
    text(topo.x,topo.y,ccep(run).ch)

    ax = gca;
    xlim([min(topo.x)-2, max(topo.x)+2])
    ylim([min(topo.y)-2, max(topo.y)+2])
    title(sprintf('CCEP responses after stimulating %s-%s',ccep(run).ch{stimnum(1)},ccep(run).ch{stimnum(2)}))

    ax.YDir = 'reverse';
    ax.YTick = [];
    ax.XTick = [];
    ax.XColor = 'none';
    ax.YColor = 'none';
    ax.Units = 'normalized';
    ax.Position = [0.1 0.1 0.8 0.8];
    outlabel=sprintf('Stimpair%s-%s.jpg',...
        ccep(run).ch{stimnum(1)},ccep(run).ch{stimnum(2)});

    path = fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label,runs{run});
    if ~exist([path,'/figures/'], 'dir')
        mkdir([path,'/figures/']);
    end

    saveas(gcf,[path,'/figures/',outlabel],'jpg')

end
