function visualise_gridstructure(myDataPath, ccep, agreement_parameter)

subj = [extractBetween(ccep.dataName,'sub-','/ses')];


if exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj,'_elektroden.xls']),'Sheet','matlabsjabloon');
end

% localize electrodes in grid
% Same for 10 and 2 stims
x = NaN(size(ccep.ch)); 
y = NaN(size(ccep.ch));
elecmat = NaN(size(elec));
topo=struct;

for i=1:size(elec,1)
    for j=1:size(elec,2)
        if ~ismissing(elec{i,j})
            letter = regexp(elec{i,j},'[a-z,A-Z]');
            number = regexp(elec{i,j},'[1-9]');
            test1 = elec{i,j}([letter,number:end]);
            test2 = [elec{i,j}(letter),'0',elec{i,j}(number:end)];
            if sum(strcmp(ccep.ch,test1))==1
                elecmat(i,j) = find(strcmp(ccep(1).ch,test1));
                y(strcmp(ccep.ch,test1),1) = i;
                x(strcmp(ccep.ch,test1),1)= j;
            elseif sum(strcmp(ccep.ch,test2))==1
                elecmat(i,j) = find(strcmp(ccep.ch,test2));
                y(strcmp(ccep.ch,test2),1) = i;
                x(strcmp(ccep.ch,test2),1)= j;
            else
                error('Electrode is not found')
            end
        end
    end
end

topo.x = x;
topo.y = y;


%% CCEP responses to specific stimulus
%%% dit is nog het script om locatie van de 1 en 2'en te laten zien.
%%% dus nog niet het script om de in degree of outdegree te laten zien.
    close all;
         
    % For all stims
    figure1 = figure('Position',[508,122,889,1052]);
    axes1 = axes('Parent',figure1,'Position',[0.1 0.5 0.8 0.4]);
    hold(axes1,'on');
    plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
    xlim([min(topo.x)-2, max(topo.x)+2])
    ylim([min(topo.y)-2, max(topo.y)+2])
    axes1.YDir = 'reverse';
    axes1.YTick = [];
    axes1.XTick = [];
    axes1.XColor = 'none';
    axes1.YColor = 'none';
    axes1.Units = 'normalized';
   
    for i = 1:length(ccep.ch)
     % For all stims
        if ismember(i, agreement_parameter.highest_ind_10, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
         Ind = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','m','MarkerEdgeColor','k');
       
       elseif ismember(i,agreement_parameter.highest_outd_10,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
           Outd = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','g','MarkerEdgeColor','k');
       end
       
       if ismember(i, agreement_parameter.highest_ind_10, 'rows') && ismember(i,agreement_parameter.highest_outd_10,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
          Both =  plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','g','MarkerEdgeColor','m','LineWidth',3);
       end  
       
    end
    
    hold(axes1,'off')  
    str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)

    title('\rm Highest indegree and outdegree scoring electrodes, all stims')
    text(((topo.x)+0.2),topo.y,ccep.ch)
    legend([Ind, Outd, Both],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_10))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_10))),'Both'},'Position',[0.005,0.5,0.15,0.12])

    
    % Plot for 2 stims      
    axes2 = axes('Parent',figure1,'Position',[0.1 0.07 0.8 0.4]);
    hold(axes2,'on');
    plot(topo.x,topo.y,'ok','Parent',axes2,'MarkerSize',15);
    xlim([min(topo.x)-2, max(topo.x)+2])
    ylim([min(topo.y)-2, max(topo.y)+2])
    axes2.YDir = 'reverse';
    axes2.YTick = [];
    axes2.XTick = [];
    axes2.XColor = 'none';
    axes2.YColor = 'none';
    axes2.Units = 'normalized';

    for i = 1:length(ccep.ch)
       if ismember(i, agreement_parameter.highest_ind_2, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
         Ind_2 = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','m','MarkerEdgeColor','k');
       
       elseif ismember(i,agreement_parameter.highest_outd_2,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
           Outd_2 = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','g','MarkerEdgeColor','k');
       end
       
       if ismember(i, agreement_parameter.highest_ind_2, 'rows') && ismember(i,agreement_parameter.highest_outd_2,'rows')               % When the electrode is highest ranked in the indegree and in outdegree
          Both_2 =  plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
               'MarkerFaceColor','g','MarkerEdgeColor','m','LineWidth',3);
       end     
    end 
        
    %hold(axes2,'off')  
    title('\rm Highest indegree and outdegree scoring electrodes, 2 stims')
    text(((topo.x)+0.2),topo.y,ccep.ch)
    outlabel=sprintf('sub-%s.jpg',...
    subj{1});
    path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/',subj{1}),'/'];
    legend([Ind_2, Outd_2, Both_2],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_2))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_2))),'Both'},'Position',[0.005,0.005,0.15,0.12])
    
    % Save figure
    if ~exist(path, 'dir')
       mkdir(path);
    end    
    saveas(gcf,[path,outlabel],'jpg')
  
    
    %% Plot lines in between electrodes to show the number of ERs evoked by that stimulation pair
    % Divide number of ERs per stimpair in 7 equal bin
    % For all stimpairs
    edges_10 = linspace(min(agreement_parameter.ERs_stimp10), max(agreement_parameter.ERs_stimp10),7);            % create 7 equally-spaced bins based on their indegree score
    bins_10 = discretize(agreement_parameter.ERs_stimp10, edges_10);
    high_bins_10 = max(bins_10);
    highest_stimp_10 = find(bins_10 == high_bins_10);
   
    % For 2 stimpairs
    edges_2 = linspace(min(agreement_parameter.ERs_stimp2), max(agreement_parameter.ERs_stimp2),7);            % create 7 equally-spaced bins based on their indegree score
    bins_2 = discretize(agreement_parameter.ERs_stimp2, edges_2);
    high_bins_2 = max(bins_2);
    highest_stimp_2 = find(bins_2 == high_bins_2);
   
   
    %axes2 = axes('Parent',figure1,'Position',[0.1 0.07 0.8 0.4]);
    %hold(axes2,'on');
   % plot(topo.x,topo.y,'ok','Parent',axes2,'MarkerSize',15);
%     xlim([min(topo.x)-2, max(topo.x)+2])
%     ylim([min(topo.y)-2, max(topo.y)+2])
%     axes2.YDir = 'reverse';
%     axes2.YTick = [];
%     axes2.XTick = [];
%     axes2.XColor = 'none';
%     axes2.YColor = 'none';
%     axes2.Units = 'normalized';

    for i = 1:length(ccep.stimchans_avg)
       if ismember(i, highest_stimp_2, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
         elec1 = ccep.stimsets_avg(i,1);
         elec2 = ccep.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'r');
         set(fig,{'LineWidth'},{3})
       % elseif onderdeel van een lagere groep, dan dus linedikte iets
       % smaller zetten.
       end     
    end 
        
    hold(axes2,'off')  
    title('\rm Highest indegree and outdegree scoring electrodes, 2 stims')
    text(((topo.x)+0.2),topo.y,ccep.ch)
    outlabel=sprintf('sub-%s.jpg',...
    subj{1});
    path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/',subj{1}),'/'];
    legend([Ind_2, Outd_2, Both_2],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_2))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_2))),'Both'},'Position',[0.005,0.005,0.15,0.12])
    
   
   
   % mogelijk ook nog iets van lijnen ertussen plotten om zo'n connectie te
   % laten zien?
   % for 1: aantal stimulatieparen
   % order aan de hand van het aantal ERs per stimpaar
   % vervolgens weer 7 van die equal Bins maken
   % stimparen uit hogere Bins dikkere line geven --> LineWitdh
   % (mogelijk dit? --> plot(topo.x,topo.y,'k--'))
   %


  
end
