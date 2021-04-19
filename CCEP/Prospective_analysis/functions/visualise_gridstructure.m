function visualise_gridstructure(myDataPath, ccep_clin, agreement_parameter)
% Display the values of the network parameters/characterisitcs on the
% patients grid-structure from the matlabsjabloon sheet in an excel file

subj = extractBetween(ccep_clin.dataName,'sub-','/ses');

if exist(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xlsx']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xlsx']),'Sheet','matlabsjabloon');
elseif exist(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xls']),'file')
    elec = readcell(fullfile(myDataPath.elec_input,[subj{1},'_ses-1_elektroden.xls']),'Sheet','matlabsjabloon');
end

% localize electrodes in grid
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
mode = {'Indegree & ERs per stimpair, Clinical SPES','Indegree & ERs per stimpair, Propofol SPES'};
figure1 = figure('Name',subj{:},'Position',[284,4,1309,1052]);
    
 for J = 1:size(mode,2)
    
    if strcmp(mode{J},'Indegree & ERs per stimpair, Clinical SPES')
        Ind = (agreement_parameter.indegreeN_Clin)';
        ERs = agreement_parameter.ERs_stimpClin;
        axes1 = axes('Parent',figure1,'Position',[0.04,0.5,0.9,0.4]);
        
    elseif strcmp(mode{J},'Indegree & ERs per stimpair, Propofol SPES')
        Ind = (agreement_parameter.indegreeN_Prop)';
        ERs = agreement_parameter.ERs_stimpProp;
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
    
    for i = 1:length(ccep_clin.stimchans_avg)
        elec1 = ccep_clin.stimsets_avg(i,1);                         % Find the electrodes in the stimpair
        elec2 = ccep_clin.stimsets_avg(i,2);
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
    
%      % Create the same colormap limits based on the highest protocol with the highest indegree
%     if max(agreement_parameter.indegreeN_Clin)>max(agreement_parameter.indegreeN_Prop)
%         caxis([0 max(agreement_parameter.indegreeN_Clin)]);
%     else
%         caxis([0 max(agreement_parameter.indegreeN_Prop)]);
%     end
%     
    hold(axes1,'on')
    str_main = sprintf('sub-%s', subj{1});
    sgtitle(str_main)
    text(((topo.x)+0.2),topo.y,ccep_clin.ch, 'FontSize',8,'FontWeight','bold')
    
 end    
    
 % Save figure
outlabel=sprintf('sub-%s_indegree_ERstimp.png',...
    subj{1});
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Visualise Gridstructure/');
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'png')

%% Plot the outdegree and the BC

mode = {'Outdegree','BC','Indegree'};

for J = 1:size(mode,2) 
    
    figure2 = figure('Name',subj{:},'Position',[284,4,1309,1052]);
    
    if strcmp(mode{J},'Outdegree')
        parclin = (agreement_parameter.outdegreeN_Clin)';
        parprop = (agreement_parameter.outdegreeN_Prop)';
        
    elseif strcmp(mode{J},'BC')
        parclin = (agreement_parameter.BCN_Clin)';
        parprop = (agreement_parameter.BCN_Prop)';
        
    elseif strcmp(mode{J},'Indegree')
        parclin = (agreement_parameter.indegreeN_Clin)';
        parprop = (agreement_parameter.indegreeN_Prop)';
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
    title(sprintf('%s, Clinical SPES',mode{J}))
    
    % Plot the parameter as a colorscale
    scatter(topo.x, topo.y, 260, parclin,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    colorbar();
    text(((topo.x)+0.2),topo.y,ccep_clin.ch, 'FontSize',8,'FontWeight','bold')
    
%     % Create the same colormap limits based on the highest protocol with
%     % the highest parameter value
%     if max(par10)>max(par2)
%         caxis([0 max(par10)]);
%     else
%         caxis([0 max(par2)]);
%     end
    
    % Prop stims
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
    title(sprintf('%s, Propofol SPES',mode{J}))
    
    % Plot the outdegree as a colorscale, combine this with the number of
    % ERS evoked per stimulation pair.
    scatter(topo.x, topo.y, 260, parprop,'filled','MarkerEdgeColor','k')
    c = hot;
    c = flipud(c);
    colormap(c);
    colorbar();
    text(((topo.x)+0.2),topo.y,ccep_clin.ch, 'FontSize',8,'FontWeight','bold')
    
    % Create the same colormap limits based on the highest protocol with
    % the highest parameter value
%     if max(par10)>max(par2)
%         caxis([0 max(par10)]);
%     else
%         caxis([0 max(par2)]);
%     end
    
    % Save figure
    outlabel=sprintf('sub-%s_%s.png',...
        subj{1},mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Visualise Gridstructure/');
    
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'png')
    
end

  close all;

%% ER's responses to specific stimulus
plot_fig = input('Do you want plot figures with all ERs per stimulation pair? [y/n] ','s');

if strcmp(plot_fig,'y')

    for stimp = 1:size(ccep_clin.stimsets_avg)                   % Number of stimulation pairs (columns)
        stimnum = ccep_clin.stimsets_avg(stimp,:);            % Stimulation pair numbers for column number (stimp)

        
        figure3= figure('Position',[284,4,1309,1052]);
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

        % for CLinicalA SPES stims
        % plot stimulation pair in black
        for chan=1:2
            indicator(1,:) = plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor','k','MarkerEdgeColor','k');
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');
        
        
        % Plot ERs in grey
        for elek = 1:length(ccep.ch)
            if ~isnan(ccep.n1_peak_sample(elek,stimp))
                indicator(2,:) =  plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor', [0.75 0.75 0.75],'MarkerEdgeColor','k');
            end
            
        end
        
        % In blue the ERs in Prop and not in clin
        for elek = 1:length(ccep.ch)
            if isnan(ccep.n1_peak_sample(elek,stimp)) && ~isnan(ccep2.n1_peak_sample(elek,stimp))
                indicator(3,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[0 0.7 1],'MarkerEdgeColor','k');
            end
        end
        
        % In magenta the ER in CLin and not in prop
        for elek = 1:length(ccep.ch)
            if isnan(ccep2.n1_peak_sample(elek,stimp)) && ~isnan(ccep.n1_peak_sample(elek,stimp))
                indicator(4,:) = plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[1 0.7 1],'MarkerEdgeColor','k');
            end
        end
        
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        title('\rm ERs responses to specific stimulus, all stims')
        
        % Save figure
        outlabel=sprintf('sub-%s, stimp %s.png',subj{1},ccep_clin.stimpnames_avg{stimp});
        path = fullfile(myDataPath.CCEPpath,'Grid_diffRes/',subj{1});
        
        if ~exist(path, 'dir')
            mkdir(path);
        end
        saveas(gcf,[path,outlabel],'png')
        
    end 
        
        
end
close all

end
