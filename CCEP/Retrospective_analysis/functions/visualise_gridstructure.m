function visualise_gridstructure(myDataPath, ccep, ccep2, agreement_parameter,plot_fig)

subj = [extractBetween(ccep.dataName,'sub-','/ses')];
scale_10 = (agreement_parameter.indegreeN_10)';
scale_2 = (agreement_parameter.indegreeN_2)';

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
% For all stimpairs
figure1 = figure('Name',subj{:},'Position',[284,4,1309,1052]);
axes1 = axes('Parent',figure1,'Position',[0.04,0.5,0.9,0.4]);
hold(axes1,'on');
plot(topo.x,topo.y,'ok','Parent',axes1,'MarkerSize',15);
xlim([min(topo.x)-1, max(topo.x)+1])
ylim([min(topo.y)-2, max(topo.y)+2])
axes1.YDir = 'reverse';                                         % Flip figure
axes1.YTick = [];                                               % Remove numbers on y-axis
axes1.XTick = [];                                               % Remove numbers on x-axis
axes1.XColor = 'none';                                          % Remove line indicating the x-axis
axes1.YColor = 'none';                                          % Remove line indicating the y-axis

edges_10 = linspace(min(agreement_parameter.ERs_stimp10), max(agreement_parameter.ERs_stimp10),7);            % create 7 equally-spaced bins based on their indegree score
bins_10 = discretize(agreement_parameter.ERs_stimp10, edges_10);
high_bins_10 = max(bins_10);
highest_10_1 = find(bins_10 == high_bins_10);
highest_10_2 = find(bins_10 == high_bins_10-1);
highest_10_3 = find(bins_10 == high_bins_10-2);

% Draw lines between the stimulation pairs with the highest number of ERs
for i = 1:length(ccep.stimchans_avg)
    if ismember(i, highest_10_1, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{5})
        fig.Color(4) = 0.4;
        
    elseif ismember(i, (highest_10_2), 'rows')
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{3})
        fig.Color(4) = 0.4;
        
    elseif ismember(i, (highest_10_3), 'rows')
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{1})
        fig.Color(4) = 0.4;
        
    end
end

% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.
scatter(topo.x, topo.y, 260, scale_10,'filled','MarkerEdgeColor','k')
c = hot;
c = flipud(c);
colormap(c);
cb = colorbar();

% Create the same colormap limits based on the highest protocol with the highest indegree
if max(scale_10)>max(scale_2)                                   
    caxis([0 max(scale_10)]);
else
    caxis([0 max(scale_2)]);
end


hold(axes1,'off')
str_main = sprintf('sub-%s', subj{1});
sgtitle(str_main)

title({'\rm Highest indegree scoring electrodes are darker green,'...
    'broader lines indicate more ERs evoked per stimulation pair, all stims'})
text(((topo.x)+0.2),topo.y,ccep.ch, 'FontSize',8)

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
for i = 1:length(ccep.stimchans_avg)
    if ismember(i, highest_2_1, 'rows')              % When the electrode is highest ranked in the indegree and in outdegree
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{5})
        fig.Color(4) = 0.4;
        
    elseif ismember(i, (highest_2_2), 'rows')
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{3})
        fig.Color(4) = 0.4;
        
    elseif ismember(i, (highest_2_3), 'rows')
        elec1 = ccep.stimsets_avg(i,1);
        elec2 = ccep.stimsets_avg(i,2);
        fig = plot([topo.x(elec1), topo.x(elec2)], [topo.y(elec2),topo.y(elec2)], 'k');
        set(fig,{'LineWidth'},{1})
        fig.Color(4) = 0.4;
        
    end
end


% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.
scatter(topo.x, topo.y, 260, scale_2,'filled','MarkerEdgeColor','k')
c = hot;
c = flipud(c);
colormap(c);
caxis([0 max(scale_2)]);
cb = colorbar();

title({'\rm Highest indegree scoring electrodes are darker green,'...
    'broader lines indicate more ERs evoked per stimulation pair, 2 stims'})
text(((topo.x)+0.2),topo.y,ccep.ch,'FontSize',8)


% Create the same colormap limits based on the highest protocol with the highest indegree
if max(scale_10)>max(scale_2)                                   
    caxis([0 max(scale_10)]);
else
    caxis([0 max(scale_2)]);
end


% Save figure
outlabel=sprintf('sub-%s_indegree_ERstimp.jpg',...
    subj{1});
path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/')];
if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')

%% Outdegree of electrodes, for all stims
% All stims
scale_2_out = (agreement_parameter.outdegreeN_2)';
scale_10_out = (agreement_parameter.outdegreeN_10)';

figure2 = figure('Name',subj{:},'Position',[284,4,1309,1052]);
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
title('\rm Highest outdegree scoring electrodes are darker, all stims')

% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.
scatter(topo.x, topo.y, 260, scale_10_out,'filled','MarkerEdgeColor','k')
c = hot;
c = flipud(c);
colormap(c);   
cb = colorbar();
text(((topo.x)+0.2),topo.y,ccep.ch,'bold')

% Create the same colormap limits based on the highest protocol with the highest outdegree
if max(scale_10_out)>max(scale_2_out)                                   
    caxis([0 max(scale_10_out)]);
else
    caxis([0 max(scale_2_out)]);
end


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
str_main = sprintf('sub-%s', subj{1});
sgtitle(str_main)
title('\rm Highest outdegree scoring electrodes are darker, 2 stims')

% Plot the indegree as a colorscale, combine this with the number of
% ERS evoked per stimulation pair.
scatter(topo.x, topo.y, 260, scale_2_out,'filled','MarkerEdgeColor','k')
c = hot;
c = flipud(c);
colormap(c);
cb = colorbar();
text(((topo.x)+0.2),topo.y,ccep.ch,'bold')

% Create the same colormap limits based on the highest protocol with the highest outdegree
if max(scale_10_out)>max(scale_2_out)                                   
    caxis([0 max(scale_10_out)]);
else
    caxis([0 max(scale_2_out)]);
end

% Save figure
outlabel=sprintf('sub-%s_outdegree.jpg',...
    subj{1});
path = [fullfile(myDataPath.CCEPpath,'Visualise_agreement/')];

if ~exist(path, 'dir')
    mkdir(path);
end
saveas(gcf,[path,outlabel],'jpg')

%% ER's responses to specific stimulus
if ~exist('plot_fig')
    plot_fig = input('Do you want plot figures with all ERs per stimulation pair? [y/n] ','s');
end

if strcmp(plot_fig,'y')
    
    for stimp = 1:size(ccep(1).stimsets_avg)                   % Number of stimulation pairs (columns)
        stimnum = ccep(1).stimsets_avg(stimp,:);               % Stimulation pair numbers for column number (stimp)
        
        % for 10 stims
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
        str_main = sprintf('sub-%s', subj{1});
        sgtitle(str_main)
        title('\rm ERs responses to specific stimulus, all stims')
        
        % plot stimulation pair in black
        for chan=1:2
            plot(topo.x(stimnum(chan)),topo.y(stimnum(chan)),'o','MarkerSize',15,...
                'MarkerFaceColor','k','MarkerEdgeColor','k')
        end
        plot([topo.x(stimnum(1)), topo.x(stimnum(2))], [topo.y(stimnum(1)),topo.y(stimnum(2))], 'k');
        
        % plot ERs in grey
        for elek = 1:length(ccep.ch)
            if ~isnan(ccep.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor', [0.75 0.75 0.75],'MarkerEdgeColor','k')
            end
            
        end
        
        % plot ERs in grey
        for elek = 1:length(ccep.ch)
            if ~isnan(ccep2.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor', [0.75 0.75 0.75],'MarkerEdgeColor','k')
            end
        end
        
        % In blue the ERs in 2 and not in 10
        for elek = 1:length(ccep.ch)
            if isnan(ccep.n1_peak_sample(elek,stimp)) && ~isnan(ccep2.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[0 0.7 1],'MarkerEdgeColor','k')
            end
        end
        
        % In magenta the ER in 10 and not in 2
        for elek = 1:length(ccep.ch)
            if isnan(ccep2.n1_peak_sample(elek,stimp)) && ~isnan(ccep.n1_peak_sample(elek,stimp))
                plot(topo.x(elek),topo.y(elek),'o','MarkerSize',15,...
                    'MarkerFaceColor',[1 0.7 1],'MarkerEdgeColor','k')
            end
        end
    end
end
end


