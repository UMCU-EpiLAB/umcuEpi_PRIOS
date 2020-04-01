%% display connections

%% pre-allocation
clear
config_CCEP

%% load all CCEP set

cfg.CCEPpath = fullfile('/Fridge/users/sifra/derivatives/CCEP/',cfg.sub_labels{:},cfg.ses_label);

files = dir(cfg.CCEPpath);
n=1; runs = cell(1);

for i=1:size(files,1)
   if contains(files(i).name ,'run-') && n==1
       loadfile = load(fullfile(cfg.CCEPpath,files(i).name,[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',files(i).name,'_CCEP.mat']));
       ccep = loadfile.ccep;
       runs{n} = files(i).name;
       n=n+1;
   elseif contains(files(i).name ,'run-') && n>1
       loadfile = load(fullfile(cfg.CCEPpath,files(i).name,[cfg.sub_labels{:},'_',cfg.ses_label,'_task-SPESclin_',files(i).name,'_CCEP.mat']));
       ccep(n) = loadfile.ccep;
       runs{n} = files(i).name;
       n=n+1;      
   end
end
 
%% load electrodes positions (xlsx/electrodes.tsv)
%cfg.proj_dirinput = '/home/dorien/Desktop/bulkstorage/db/respect-leijten/Electrodes/';
cfg.proj_dirinput = '/home/sifra/Desktop/db/Electrodes/';

subj = cfg.sub_labels{1}(5:end);

% if exist(fullfile(cfg.proj_dirinput,[subj,'_',cfg.ses_label,'_elektroden.xlsx']),'file')
%     elec = readcell(fullfile(cfg.proj_dirinput,[subj,'_',cfg.ses_label,'_elektroden.xlsx']),'Sheet','matlabsjabloon','Range',[1 1 100 100]);
% elseif exist(fullfile(cfg.proj_dirinput,[subj,'_',cfg.ses_label,'_elektroden.xls']),'file')
%     elec = readcell(fullfile(cfg.proj_dirinput,[subj,'_',cfg.ses_label,'_elektroden.xls']),'Sheet','matlabsjabloon','Range',[1 1 100 100]);
% end
    
if exist(fullfile(cfg.proj_dirinput,[subj,'_elektroden.xlsx']),'file')
    elec = readcell(fullfile(cfg.proj_dirinput,[subj,'_elektroden.xlsx']),'Sheet','matlabsjabloon','Range',[1 1 100 100]);
elseif exist(fullfile(cfg.proj_dirinput,[subj,'_elektroden.xls']),'file')
    elec = readcell(fullfile(cfg.proj_dirinput,[subj,'_elektroden.xls']),'Sheet','matlabsjabloon','Range',[1 1 100 100]);
end

% localize electrodes in grid
x = NaN(size(ccep(1).ch)); 
y = NaN(size(ccep(1).ch));
elecmat = NaN(size(elec));
topo=struct;
%%
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
%%
% for i=1:size(elec,1)
%     for j=1:size(elec,2)
%         if ~ismissing(elec{i,j})
%                     
%             letter = regexp(string(elec{i,j}),'[a-z,A-Z]');
%             number = regexp(string(elec{i,j}),'[1-9]');
%             
%             test1 = elec{i,j}([letter,number:end]);
%             test2 = [elec{i,j}(letter),'0',elec{i,j}(number:end)];
%             
%             if sum(strcmp(ccep(1).ch,test1))==1
%                 elecmat(i,j) = find(strcmp(ccep(1).ch,test1));
%                 y(strcmp(ccep(1).ch,test1),1) = i;
%                 x(strcmp(ccep(1).ch,test1),1)= j;
%             elseif sum(strcmp(ccep(1).ch,test2))==1
%                 elecmat(i,j) = find(strcmp(ccep(1).ch,test2));
%                 y(strcmp(ccep(1).ch,test2),1) = i; 
%                 x(strcmp(ccep(1).ch,test2),1)= j;
%             else
%                 error('Electrode is not found')
%             end
%         end
%     end
% end
%%

topo.x =x;
topo.y=y;


%% CCEP responses to specific stimulus

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
    
    path = fullfile(cfg.CCEPpath,runs{run});
    if ~exist([path,'/figures/'], 'dir')
        mkdir([path,'/figures/']);
    end
    
    saveas(gcf,[path,'/figures/',outlabel],'jpg')
    
end
end