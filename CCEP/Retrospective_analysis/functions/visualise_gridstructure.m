function visualise_gridstructure(myDataPath, ccep, ccep2, agreement_parameter,plot_fig)

subj = extractBetween(ccep.dataName,'sub-','/ses');
% scale_10 = (agreement_parameter.indegreeN_10)';
% scale_2 = (agreement_parameter.indegreeN_2)';

if exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_elektroden.xls']),'Sheet','matlabsjabloon');
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


%% Indegree of electrodes and ERs per stimulation pair, for all stims

mode = {'Indegree & ERs per stimpair, all stimuli','Indegree & ERs per stimpair, 2 stimuli'};
figure1 = figure('Name',subj{:},'Position',[284,4,1309,1052]);

for J = 1:size(mode,2)
    
    if strcmp(mode{J},'Indegree & ERs per stimpair, all stimuli')
        Ind = (agreement_parameter.indegreeN_10)';
        ERs = agreement_parameter.ERs_stimp10;
        axes1 = axes('Parent',figure1,'Position',[0.04,0.5,0.9,0.4]);
        
    elseif strcmp(mode{J},'Indegree & ERs per stimpair, 2 stimuli')
        Ind = (agreement_parameter.indegreeN_2)';
        ERs = agreement_parameter.ERs_stimp2;
        axes1 = axes('Parent',figure1,'Position',[0.04,0.07,0.9,0.4]);
        
    end
    
    hold(axes1,'on');
    plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
    xlim([min(topo.x)-1, max(topo.x)+1])
    ylim([min(topo.y)-0.5, max(topo.y)+0.5])
    axes1.YDir = 'reverse';                                         % Flip figure
    axes1.YTick = [];                                               % Remove numbers on y-axis
    axes1.XTick = [];                                               % Remove numbers on x-axis
    axes1.XColor = 'none';                                          % Remove line indicating the x-axis
    axes1.YColor = 'none';                                          % Remove line indicating the y-axis
    title(sprintf('%s',mode{J}))
    
    % Draw lines between stimpairs, thickness of the line indicates the number
    % of ERs per stimpair. Thicker = more ERs
    lineWidth_set = linspace(min(ERs), max(ERs),6);                 % devide the number of ERs to fit the lineWidth
    bins_lineWidth = discretize(ERs,lineWidth_set);
    
    for i = 1:length(ccep.stimchans_avg)
        elec1 = ccep.stimsets_avg(i,1);                         % Find the electrodes in the stimpair
        elec2 = ccep.stimsets_avg(i,2);
        if ERs(i) > 0                                           % Stimpairs with 0 ERs do not get a line.
            
            % horizontal lines between electrodes
            line_stim_x = line([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'Color',[0.6 0.6 0.6]-0.05*bins_lineWidth(i),'LineWidth',bins_lineWidth(i));
            line_stim_x.Color(4) = 0.6;              % Transparency
            
            % Vertical lines between electrodes
            line_stim_y = line([topo.x(elec1), topo.x(elec1)], [topo.y(elec1),topo.y(elec2)],'Color',[0.6 0.6 0.6]-0.05*bins_lineWidth(i),'LineWidth',bins_lineWidth(i));
            line_stim_y.Color(4) = 0.6;
        end
        
    end
    
    % Plot the indegree as a colorscale, combine this with the number of
    % ERS evoked per stimulation pair.
    scatter(topo.x, topo.y, 260, Ind,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    colorbar();
    
    % Create the same colormap limits based on the highest protocol with the highest indegree
%     if max(agreement_parameter.indegreeN_10)>max(agreement_parameter.indegreeN_2)
%         caxis([0 max(agreement_parameter.indegreeN_10)]);
%     else
%         caxis([0 max(agreement_parameter.indegreeN_2)]);
%     end
%     
    
    hold(axes1,'on')
    str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)
    text(((topo.x)+0.2),topo.y,ccep.ch, 'FontSize',8,'FontWeight','bold')
    
end

% Save figure
outlabel=sprintf('sub-%s_indegree_ERstimp.jpg',...
    subj{1});
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')



%% Plot the outdegree and the BC

mode = {'Outdegree','BC','Indegree'};

for J = 1:size(mode,2)
    
    figure2 = figure('Name',subj{:},'Position',[284,4,1309,1052]);
    
    if strcmp(mode{J},'Outdegree')
        par10 = (agreement_parameter.outdegreeN_10)';
        par2 = (agreement_parameter.outdegreeN_2)';
        
    elseif strcmp(mode{J},'BC')
        par10 = (agreement_parameter.BCN_10)';
        par2 = (agreement_parameter.BCN_2)';
   
    elseif strcmp(mode{J},'Indegree')
        par10 = (agreement_parameter.indegreeN_10)';
        par2 = (agreement_parameter.indegreeN_2)';
        
    end
    
    axes3 = axes('Parent',figure2,'Position',[0.04,0.5,0.9,0.4]);
    hold(axes3,'on');
    plot(topo.x,topo.y,'ok','Parent',axes3,'MarkerSize',15);
    xlim([min(topo.x)-1, max(topo.x)+1])
    ylim([min(topo.y)-0.5, max(topo.y)+0.5])
    axes3.YDir = 'reverse';
    axes3.YTick = [];
    axes3.XTick = [];
    axes3.XColor = 'none';
    axes3.YColor = 'none';
    title(sprintf('%s, all stimuli',mode{J}))
    
    % Plot the parameter as a colorscale
    scatter(topo.x, topo.y, 260, par10,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    colorbar();
    text(((topo.x)+0.2),topo.y,ccep.ch, 'FontSize',8,'FontWeight','bold')
    
%     % Create the same colormap limits based on the highest protocol with
%     % the highest parameter value
%     if max(par10)>max(par2)
%         caxis([0 max(par10)]);
%     else
%         caxis([0 max(par2)]);
%     end
    
    
    % 2 stims
    axes4 = axes('Parent',figure2,'Position',[0.04,0.07,0.9,0.4]);
    hold(axes4,'on');
    plot(topo.x,topo.y,'ok','Parent',axes4,'MarkerSize',15);
    xlim([min(topo.x)-1, max(topo.x)+1])
    ylim([min(topo.y)-0.5, max(topo.y)+0.5])
    axes4.YDir = 'reverse';
    axes4.YTick = [];
    axes4.XTick = [];
    axes4.XColor = 'none';
    axes4.YColor = 'none';
    str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)
    title(sprintf('%s, 2 stimuli',mode{J}))
    
    % Plot the outdegree as a colorscale, combine this with the number of
    % ERS evoked per stimulation pair.
    scatter(topo.x, topo.y, 260, par2,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    colorbar();
    text(((topo.x)+0.2),topo.y,ccep.ch, 'FontSize',8,'FontWeight','bold')
    
    % Create the same colormap limits based on the highest protocol with
    % the highest parameter value
%     if max(par10)>max(par2)
%         caxis([0 max(par10)]);
%     else
%         caxis([0 max(par2)]);
%     end
    
    % Save figure
    outlabel=sprintf('sub-%s_%s.jpg',...
        subj{1},mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/');
    
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')
    
end



%% ER's responses to specific stimulus
if ~exist('plot_fig','var')
    plot_fig = input('Do you want plot figures with all ERs per stimulation pair & SOZ? [y/n] ','s');
end

if strcmp(plot_fig,'y')
    
    for stimp = 1:size(ccep(1).stimsets_avg)                   % Number of stimulation pairs (columns)
        stimnum = ccep(1).stimsets_avg(stimp,:);               % Stimulation pair numbers for column number (stimp)
        
        figure3 = figure('Position',[280,400,1300,500]);
        axes1 = axes('Parent',figure3, 'Position', [0.04,0.014,0.9,0.886]);
        hold(axes1,'on');
        plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
        xlim([min(topo.x)-1, max(topo.x)+1])
        ylim([min(topo.y)-2, max(topo.y)+2])
        axes1.YDir = 'reverse';
        axes1.YTick = [];
        axes1.XTick = [];
        axes1.XColor = 'none';
        axes1.YColor = 'none';
        text(((topo.x)+0.2),topo.y,ccep.ch,'bold')
        
        
        % for 10 stims
        % plot stimulation pair in black
        for chan=1:2
            indicator(1,:) = plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor','k','MarkerEdgeColor','k');
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');
        
        % plot ERs in grey
        for elek = 1:length(ccep.ch)
            if ~isnan(ccep.n1_peak_sample(elek,stimp))
                indicator(2,:) =  plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor', [0.75 0.75 0.75],'MarkerEdgeColor','k');
            end
            
        end
        
        % In blue the ERs in 2 and not in 10
        for elek = 1:length(ccep.ch)
            if isnan(ccep.n1_peak_sample(elek,stimp)) && ~isnan(ccep2.n1_peak_sample(elek,stimp))
                indicator(3,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[0 0.7 1],'MarkerEdgeColor','k');
            end
        end
        
        % In magenta the ER in 10 and not in 2
        for elek = 1:length(ccep.ch)
            if isnan(ccep2.n1_peak_sample(elek,stimp)) && ~isnan(ccep.n1_peak_sample(elek,stimp))
                indicator(4,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[1 0.7 1],'MarkerEdgeColor','k');
            end
        end
        
        %         legend(indicator,'Stimpair', 'ER in both','ER in 2','ER in 10')
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        title('\rm ERs responses to specific stimulus, all stims')
        
        % Save figure
        outlabel=sprintf('sub-%s, stimp %s.jpg',subj{1},ccep.stimpnames_avg{stimp});
        path = fullfile(myDataPath.CCEPpath,'Grid_diffRes/',subj{1});
        
        if ~exist(path, 'dir')
            mkdir(path);
        end
        saveas(gcf,[path,outlabel],'jpg')
        
    end
    
    
    
end
close all


%% Plot SOZ
%  close all
%
%  if ~exist('plot_fig')
%     plot_fig = input('Do you want plot figures with all ERs per stimulation pair & SOZ? [y/n] ','s');
%  end
%
% mode = {'Indegree 10','Indegree 2','Outdegree 10','Outdegree 2'};
%
% if strcmp(plot_fig,'y')
%
%     for i = 1:size(mode,2)
%         % Plot the SOZ, this is equal for all different modes.
%         figure4 = figure('Position',[284,4,1309,1052]);
%         axes5 = axes('Parent',figure4, 'Position', [0.04,0.5,0.9,0.4]);
%         hold(axes5,'on');
%         plot(topo.x,topo.y,'ok','Parent',axes5,'MarkerSize',15);
%         xlim([min(topo.x)-1, max(topo.x)+1])
%         ylim([min(topo.y)-2, max(topo.y)+2])
%         axes5.YDir = 'reverse';
%         axes5.YTick = [];
%         axes5.XTick = [];
%         axes5.XColor = 'none';
%         axes5.YColor = 'none';
%         text(((topo.x)+0.2),topo.y,ccep.ch,'bold')
%         hold on
%
%         str_main = sprintf('sub-%s', subj{1});
%         sgtitle(str_main)
%         title('\rm ELectrodes in SOZ')
%
%         % In blue electrodes in SOZ
%         for elek = 1:length(ccep.ch)
%             if ismember(ccep.SOZ(elek),'yes')
%                SOZ = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
%                     'MarkerFaceColor',[0 0.2 1],'MarkerEdgeColor','k');
%             end
%         end
%
%         legend(SOZ,'Electrode in SOZ')
%
%
%         % Per loopturn a differnt parameters is compared with the SOZ
%         if strcmp(mode{i},'Indegree 10')
%             val = scale_10;
%         elseif strcmp(mode{i},'Indegree 2')
%             val = scale_2;
%         elseif strcmp(mode{i},'Outdegree 10')
%             val = scale_10_out;
%         elseif strcmp(mode{i},'Outdegree 2')
%             val = scale_2_out;
%         end
%
%         % Plot the parameters to compare with the SOZ
%         axes6 = axes('Parent',figure4,'Position', [0.04,0.07,0.9,0.4]);
%         hold(axes6,'on');
%         plot(topo.x,topo.y,'ok','Parent',axes6,'MarkerSize',15);
%         xlim([min(topo.x)-1, max(topo.x)+1])
%         ylim([min(topo.y)-2, max(topo.y)+2])
%         axes6.YDir = 'reverse';
%         axes6.YTick = [];
%         axes6.XTick = [];
%         axes6.XColor = 'none';
%         axes6.YColor = 'none';
%         text(((topo.x)+0.2),topo.y,ccep.ch,'bold')
%         hold on
%
%         str_main = sprintf('sub-%s', subj{1});
%         sgtitle(str_main)
%         title(sprintf(' Electrodes with highest %s stimuli',mode{i}))
%
%
%         % Determine the number of SOZ electrodes
%         size_SOZ = sum(ismember(ccep.SOZ,'yes'));                                % plot the same number of electrodes of the highest scoring
%         x =maxk(val,size_SOZ);
%
%
%         % In light red electrodes electrodes with high parameter value
%         for elek = 1:length(ccep.ch)
%             if ismember(val(elek), x)
%                val_high(1,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
%                     'MarkerFaceColor',[1 0.8 0.8],'MarkerEdgeColor','k');
%             end
%         end
%
%         hold on
%
%          % In red electrodes iin both SOZ and indegree
%         for elek = 1:length(ccep.ch)
%             if ismember(val(elek), x) && ismember(ccep.SOZ(elek),'yes')
%                val_high(2,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
%                     'MarkerFaceColor','r','MarkerEdgeColor','k');
%             end
%         end
%
%         legend(val_high,sprintf('Highest %s',mode{i}),'In both')
%
%
%          % Save figure
%         outlabel=sprintf('sub-%s_%s.jpg',subj{1},mode{i});
%         path = [fullfile(myDataPath.CCEPpath,'SOZ/')];
%
%         if ~exist(path, 'dir')
%             mkdir(path);
%         end
%         saveas(gcf,[path,outlabel],'jpg')
%     end
%  end

end


