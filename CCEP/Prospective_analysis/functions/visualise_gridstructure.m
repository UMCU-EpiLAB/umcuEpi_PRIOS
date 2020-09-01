function visualise_gridstructure(myDataPath, ccep_clin, ccep_prop, agreement_parameter)

subj = [extractBetween(ccep_clin.dataName,'sub-','/ses')];


if exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xls']),'Sheet','matlabsjabloon');
end

% localize electrodes in grid
% Same for 10 and 2 stims
x = NaN(size(ccep_clin.ch)); 
y = NaN(size(ccep_clin.ch));
elecmat = NaN(size(elec));
topo=struct;

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

topo.x = x;
topo.y = y;


%% Indegree of electrodes and ERs per stimulation pair, for all stims
    close all;
         
    % For all stimpairs
    figure1 = figure('Position',[284,4,1309,1052]);
    axes1 = axes('Parent',figure1,'Position',[0.04,0.5,0.9,0.4]);
    hold(axes1,'on');
    plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
    xlim([min(topo.x)-1, max(topo.x)+1])
    ylim([min(topo.y)-2, max(topo.y)+2])
    axes1.YDir = 'reverse';
    axes1.YTick = [];
    axes1.XTick = [];
    axes1.XColor = 'none';
    axes1.YColor = 'none';
    axes1.Units = 'normalized';
    
    edges_10 = linspace(min(agreement_parameter.ERs_stimp10), max(agreement_parameter.ERs_stimp10),7);            % create 7 equally-spaced bins based on their indegree score
    bins_10 = discretize(agreement_parameter.ERs_stimp10, edges_10);
    high_bins_10 = max(bins_10);
    highest_10_1 = find(bins_10 == high_bins_10);
    highest_10_2 = find(bins_10 == high_bins_10-1);
    highest_10_3 = find(bins_10 == high_bins_10-2);

   
% Draw lines between the stimulation pairs with the highest number of ERs
    for i = 1:length(ccep_clin.stimchans_avg)
       if ismember(i, highest_10_1, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{5})
         fig.Color(4) = 0.4;

         
       elseif ismember(i, (highest_10_2), 'rows') 
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{3})
         fig.Color(4) = 0.4;
         
      elseif ismember(i, (highest_10_3), 'rows') 
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{1})
         fig.Color(4) = 0.4;
         
       end     
    end
    
    % Plot the indegree as a colorscale, combine this with the number of
    % ERS evoked per stimulation pair.
    scale_10 = (agreement_parameter.indegreeN_10)';
    cdata = scale_10;
   
    scatter(topo.x, topo.y, 260, cdata,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    cb = colorbar();
    
    hold(axes1,'off')  
    str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)

    title({'\rm Highest indegree scoring electrodes are darker green,'...
        'broader lines indicate more ERs evoked per stimulation pair, all stims'})
    text(((topo.x)+0.2),topo.y,ccep_clin.ch, 'FontSize',8)

    

    
    
%     for i = 1:length(ccep.ch)
%      % For all stims
%         if ismember(i, agreement_parameter.highest_ind_10, 'rows')              % When the electrode is highest ranked in the indegree 
%          Ind = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','m','MarkerEdgeColor','k');
%        
%        elseif ismember(i,agreement_parameter.highest_outd_10,'rows')               % When the electrode is highest ranked in the  outdegree
%            Outd = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','g','MarkerEdgeColor','k');
%            
%         elseif ~ismember(i, agreement_parameter.highest_ind_10, 'rows') && ~ismember(i,agreement_parameter.highest_outd_10,'rows')           
%            None=  plot(topo.x(i),topo.y(i),'ok','Parent',axes1,'MarkerSize',15,...
%                'MarkerFaceColor','white','MarkerEdgeColor','k');
%         end
%     end
%     
%     legend([Ind, Outd],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_10))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_10)))},'Position',[0.005,0.5,0.15,0.12])
% 
%       for i = 1:length(ccep.ch)  
%        if ismember(i, agreement_parameter.highest_ind_10, 'rows') && ismember(i,agreement_parameter.highest_outd_10,'rows')             % When the electrode is highest ranked in the indegree and in outdegree
%           Both =  plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','g','MarkerEdgeColor','m','LineWidth',3);
%           legend([Ind, Outd, Both],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_10))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_10))),'Both'},'Position',[0.005,0.5,0.15,0.12])
%        end   
%       end
str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)

    title({'\rm Highest indegree scoring electrodes are darker green,'...
        'broader lines indicate more ERs evoked per stimulation pair, all stims'})
    
   
    
%% Indegree of electrodes and ERs per stimulation pair, for 2 stims    
    axes2 = axes('Parent',figure1,'Position',[0.04,0.07,0.9,0.4]);          %[0.04,0.5,0.9,0.4]
    hold(axes2,'on');
    plot(topo.x,topo.y,'ok','Parent',axes2,'MarkerSize',15);
    xlim([min(topo.x)-1, max(topo.x)+1])
    ylim([min(topo.y)-2, max(topo.y)+2])
    axes2.YDir = 'reverse';
    axes2.YTick = [];
    axes2.XTick = [];
    axes2.XColor = 'none';
    axes2.YColor = 'none';
    axes2.Units = 'normalized';
    
    % Plot lines between electrodes to show the number of ERs evoked by that stimulation pair
    % Divide number of ERs per stimpair in 7 equal bin   
    % For 2 stimpairs
    edges_2 = linspace(min(agreement_parameter.ERs_stimp2), max(agreement_parameter.ERs_stimp2),7);            % create 7 equally-spaced bins based on their indegree score
    bins_2 = discretize(agreement_parameter.ERs_stimp2, edges_2);
    high_bins_2 = max(bins_2);
    highest_2_1= find(bins_2 == high_bins_2);
    highest_2_2 = find(bins_2 == high_bins_2-1);
    highest_2_3 = find(bins_2 == high_bins_2-2);
   
% Draw lines between the stimulation pairs with the highest number of ERs
    for i = 1:length(ccep_clin.stimchans_avg)
       if ismember(i, highest_2_1, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{5})
         fig.Color(4) = 0.4;
         
       elseif ismember(i, (highest_2_2), 'rows') 
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{3})
         fig.Color(4) = 0.4;

      elseif ismember(i, (highest_2_3), 'rows') 
         elec1 = ccep_clin.stimsets_avg(i,1);
         elec2 = ccep_clin.stimsets_avg(i,2);
         fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
         set(fig,{'LineWidth'},{1})
         fig.Color(4) = 0.4;

       end     
    end 
   
    
    % Plot the indegree as a colorscale, combine this with the number of
    % ERS evoked per stimulation pair.
    scale_2 = (agreement_parameter.indegreeN_2)';
    cdata_2 = scale_2;
   
    scatter(topo.x, topo.y, 260, cdata_2,'filled','MarkerEdgeColor','k')
    %c = [0.83 0.14 0.14;  1.00 0.54 0.00;  0.47 0.25 0.80;   0.25 0.80 0.54];
    %c = [0 0.5 1; 0.5 0 1; 0.7 0.7 0.7]
    c = hot;
    c = flipud(c);
    colormap(c);
    cb = colorbar();
    
    title({'\rm Highest indegree scoring electrodes are darker green,'...
        'broader lines indicate more ERs evoked per stimulation pair, 2 stims'})
    text(((topo.x)+0.2),topo.y,ccep_clin.ch,'FontSize',8)

    
% Save figure 
outlabel=sprintf('sub-%s.jpg',...
subj{1});
path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/')];     
outlabel=sprintf('sub-%s_indegree_ERstimp.jpg',...
subj{1});
path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/')];
if ~exist(path, 'dir')
   mkdir(path);
end    
saveas(gcf,[path,outlabel],'jpg')


%     for i = 1:length(ccep.ch)
%        if ismember(i, agreement_parameter.highest_ind_2, 'rows')                  
%            Ind_2 = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','m','MarkerEdgeColor','k');
%        
%        elseif ismember(i,agreement_parameter.highest_outd_2,'rows')               
%            Outd_2 = plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','g','MarkerEdgeColor','k');
%            
%        elseif ~ismember(i, agreement_parameter.highest_ind_2, 'rows') && ~ismember(i,agreement_parameter.highest_outd_2,'rows')           
%             None=  plot(topo.x(i),topo.y(i),'ok','Parent',axes2,'MarkerSize',15,...
%                 'MarkerFaceColor','white','MarkerEdgeColor','k');          
%        end
%     end
%     legend([Ind_2, Outd_2],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_2))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_2))),'Both'},'Position',[0.005,0.005,0.15,0.12])
% 
%      
%     for i = 1:length(ccep.ch)         
%        if ismember(i, agreement_parameter.highest_ind_2, 'rows') && ismember(i,agreement_parameter.highest_outd_2,'rows')               
%            Both_2 =  plot(topo.x(i),topo.y(i),'o','MarkerSize',15,...
%                'MarkerFaceColor','g','MarkerEdgeColor','m','LineWidth',3);
%            legend([Ind_2, Outd_2, Both_2],{(sprintf('Indegree (%d)',length(agreement_parameter.highest_ind_2))),(sprintf('Outdegree (%d)',length(agreement_parameter.highest_outd_2))),'Both'},'Position',[0.005,0.005,0.15,0.12])
%        end     
%     end 
        
    %hold(axes2,'off')  

    
%% Outdegree of electrodes, for all stims

% All stims
figure2 = figure('Position',[284,4,1309,1052]);
axes3 = axes('Parent',figure2,'Position',[0.04,0.5,0.9,0.4]);
hold(axes3,'on');
plot(topo.x,topo.y,'ok','Parent',axes3,'MarkerSize',15);
xlim([min(topo.x)-2, max(topo.x)+2])
ylim([min(topo.y)-2, max(topo.y)+2])
axes3.YDir = 'reverse';
axes3.YTick = [];
axes3.XTick = [];
axes3.XColor = 'none';
axes3.YColor = 'none';
axes3.Units = 'normalized';
title('\rm Highest outdegree scoring electrodes are darker, all stims')

% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.
scale_10_out = (agreement_parameter.outdegreeN_10)';
cdata_out = scale_10_out;
scatter(topo.x, topo.y, 260, cdata_out,'filled')
c = hot;
c = flipud(c);
colormap(c);
cb = colorbar();
text(((topo.x)+0.2),topo.y,ccep_clin.ch,'bold')


% 2 stims
axes4 = axes('Parent',figure2,'Position',[0.04,0.07,0.9,0.4]);
hold(axes4,'on');
plot(topo.x,topo.y,'ok','Parent',axes4,'MarkerSize',15);
xlim([min(topo.x)-2, max(topo.x)+2])
ylim([min(topo.y)-2, max(topo.y)+2])
axes4.YDir = 'reverse';
axes4.YTick = [];
axes4.XTick = [];
axes4.XColor = 'none';
axes4.YColor = 'none';
axes4.Units = 'normalized';
str_main = sprintf('sub-%s', subj{1});
sgtitle(str_main)
title('\rm Highest outdegree scoring electrodes are darker, 2 stims')
    
% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.         
scale_2_out = (agreement_parameter.outdegreeN_2)';
cdata_2 = scale_2_out;
scatter(topo.x, topo.y, 260, cdata_2,'filled')
c = hot;
c = flipud(c);
colormap(c);
cb = colorbar();
text(((topo.x)+0.2),topo.y,ccep_clin.ch,'bold')


% Save figure    
outlabel=sprintf('sub-%s_outdegree.jpg',...
subj{1});
path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/')];

if ~exist(path, 'dir')
   mkdir(path);
end    
saveas(gcf,[path,outlabel],'jpg')
   


%% ER's responses to specific stimulus
plot_fig = input('Do you want plot figures with all ERs per stimulation pair? [y/n] ','s');

if strcmp(plot_fig,'y')

    for stimp = 1:size(ccep_clin(1).stimsets_avg)                   % Number of stimulation pairs (columns)
        stimnum = ccep_clin(1).stimsets_avg(stimp,:);            % Stimulation pair numbers for column number (stimp)

        % for 10 stims
        figure3= figure('Position',[284,4,1309,1052]);
        %axes1 = axes('Parent',figure3,'Position', [0.05, 0.69, 0.92, 0.27]);

        axes1 = axes('Parent',figure3,'Position',[0.04,0.5,0.9,0.4]);
        hold(axes1,'on');
        plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
        xlim([min(topo.x)-1, max(topo.x)+1])
        ylim([min(topo.y)-2, max(topo.y)+2])
        axes1.YDir = 'reverse';
        axes1.YTick = [];
        axes1.XTick = [];
        axes1.XColor = 'none';
        axes1.YColor = 'none';
        axes1.Units = 'normalized';
        text(((topo.x)+0.2),topo.y,ccep_clin.ch,'bold')
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        title('\rm ERs responses to specific stimulus, all stims')   

        % plot stimulation pair in yellow
        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor',[1 0.9 0],'MarkerEdgeColor','k')
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');    


        for elek = 1:length(ccep_clin.ch)
            if ~isnan(ccep_clin.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                'MarkerFaceColor',[0 0.7 0],'MarkerEdgeColor','k')
            end

        end        


        % For 2 stims
        axes4 = axes('Parent',figure3,'Position',[0.04,0.07,0.9,0.4]);
        %axes4 = axes('Parent',figure3, 'Position', [0.05, 0.37, 0.92, 0.27]);

        hold(axes4,'on');
        plot(topo.x,topo.y,'ok','Parent',axes4,'MarkerSize',15);
        xlim([min(topo.x)-2, max(topo.x)+2])
        ylim([min(topo.y)-2, max(topo.y)+2])
        axes4.YDir = 'reverse';
        axes4.YTick = [];
        axes4.XTick = [];
        axes4.XColor = 'none';
        axes4.YColor = 'none';
        axes4.Units = 'normalized';
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        text(((topo.x)+0.2),topo.y,ccep_clin.ch,'bold')
        title('\rm  ERs responses to specific stimulus, 2 stims')

        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor',[1 0.9 0],'MarkerEdgeColor','k')
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');    


        for elek = 1:length(ccep_clin.ch)
            if ~isnan(ccep_prop.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                'MarkerFaceColor',[0 0.7 0],'MarkerEdgeColor','k')
            end
        end        



        % de elektroden die niet in beide zitten ander kleurtje geven
        figure4= figure('Position',[280,400,1300,500]);
        axes3 = axes('Parent',figure4,'Position',[0.04,0.014,0.9,0.886]);
    %    axes3 = axes('Parent',figure3, 'Position', [0.05, 0.05, 0.92, 0.27]);

        hold(axes3,'on');
        plot(topo.x,topo.y,'ok','MarkerSize',15);
        xlim([min(topo.x)-1, max(topo.x)+1])
        ylim([min(topo.y)-2, max(topo.y)+2])
        axes3.YDir = 'reverse';
        axes3.YTick = [];
        axes3.XTick = [];
        axes3.XColor = 'none';
        axes3.YColor = 'none';
        axes3.Units = 'normalized';
        text(((topo.x)+0.2),topo.y,ccep_clin.ch,'bold')
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        title('\rm ERs differently detected in the 2 stimuli or 10 stimuli protocol')  

        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor',[1 0.9 0],'MarkerEdgeColor','k')
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');    


        for elek = 1:length(ccep_clin.ch)
            if isnan(ccep_clin.n1_peak_sample(elek,stimp)) && ~isnan(ccep_prop.n1_peak_sample(elek,stimp)) 
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                'MarkerFaceColor',[0 0.7 1],'MarkerEdgeColor','k')
            end
        end 


         for elek = 1:length(ccep_clin.ch)
            if isnan(ccep_prop.n1_peak_sample(elek,stimp)) && ~isnan(ccep_clin.n1_peak_sample(elek,stimp)) 
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                'MarkerFaceColor',[1 0.7 1],'MarkerEdgeColor','k')
            end
         end 
    end
end



        %         
%         
%         % plot electrodes showing CCEPs in green (CCEP = 2 because 2 and 10 stimulations are compared)
%         chan = find(resp==2);
%         plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
%             'MarkerFaceColor','g','MarkerEdgeColor','k')
%         
%         % plot electrodes showing CCEPs in green (CCEP = 2 because 2 and 10 stimulations are compared)
%         chan = find(resp==0);
%         plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
%             'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','k')
%                 
%         % plot electrodes showing CCEPs in one of the two stimulations in red 
%         chan = find(resp==1);
%         plot(topo.x(chan),topo.y(chan),'o','MarkerSize',15,...
%             'MarkerFaceColor','r','MarkerEdgeColor','k')
%        
%         % plot stimulation pair in yellow
%         for chan=1:2
%             plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
%                 'MarkerFaceColor','y','MarkerEdgeColor','k')
%         end
%         
%         hold off
%         
%         % add electrode names
%         text(topo.x,topo.y,ccep(1).ch)
%         
%         ax = gca;
%         xlim([min(topo.x)-2, max(topo.x)+2])
%         ylim([min(topo.y)-2, max(topo.y)+2])
%         title(sprintf('CCEP responses after stimulating %s-%s', ccep(1).ch{stimnum(1)}, ccep(1).ch{stimnum(2)}))
%         
%         ax.YDir = 'reverse';
%         ax.YTick = [];
%         ax.XTick = [];
%         ax.XColor = 'none';
%         ax.YColor = 'none';
%         ax.Units = 'normalized';
%         ax.Position = [0.1 0.1 0.8 0.8];
%         outlabel=sprintf('Stimpair%s-%s.jpg',...
%             ccep(1).ch{stimnum(1)},ccep(1).ch{stimnum(2)});
%         
%         path = fullfile(myDataPath.CCEPpath,cfg.sub_labels{:},cfg.ses_label, cfg.run_label{:});
%         if ~exist([path,'/figures/'], 'dir')
%             mkdir([path,'/figures/']);
%         end
%         
%         saveas(gcf,[path,'/figures/',outlabel],'jpg')

        
        
        

