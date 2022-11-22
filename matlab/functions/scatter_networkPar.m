function  scatter_networkPar(dataBase,myDataPath)
% Visualise the agreement in a scatter plot of the absolute values
% The significance of the absolute values is determined with the Wilcoxon
% signed rank test and this does not result in an Rho value. 
% When statistically significant --> there if a difference between the 2
% measurements.

mode = {'Indegree','Outdegree','Betweenness Centrality'};
fig = figure('Position',[302,17,1110,1039]);

% Remove patients with too low interobserver agreement
skip_pat = zeros(size(dataBase,2),1);

% Remove patients that have an interobserver agreement lower than 0.6
for i = 1:size(dataBase,2)      
   if dataBase(i).ccep_clin.Ckappa < 0.6 && dataBase(i).ccep_prop.Ckappa <0.6
       skip_pat(i,:) = 1;
   end
end

dataBase(skip_pat==1) = [];

for J = 1:size(mode,2)

    clin = NaN(6,80);
    prop = NaN(6,80);

            if strcmp(mode{J},'Indegree')
               
                for b = 1:size(dataBase,2)
                  clin(b,1:size(dataBase(b).agreement_parameter.indegreeN_Clin,2)) = dataBase(b).agreement_parameter.indegreeN_Clin;
                  prop(b,1:size(dataBase(b).agreement_parameter.indegreeN_Prop,2)) = dataBase(b).agreement_parameter.indegreeN_Prop;
            
                end
            
                
            elseif strcmp(mode{J},'Outdegree')
                
                for b = 1:size(dataBase,2)
                    clin(b,1:size(dataBase(b).agreement_parameter.outdegreeN_Clin,2)) = dataBase(b).agreement_parameter.outdegreeN_Clin;
                    prop(b,1:size(dataBase(b).agreement_parameter.outdegreeN_Prop,2)) = dataBase(b).agreement_parameter.outdegreeN_Prop;            
                end
    
            elseif strcmp(mode{J},'Betweenness Centrality')
                
                for b = 1:size(dataBase,2)
                    clin(b,1:size(dataBase(b).agreement_parameter.BCN_Clin,2)) = dataBase(b).agreement_parameter.BCN_Clin;
                    prop(b,1:size(dataBase(b).agreement_parameter.BCN_Prop,2)) = dataBase(b).agreement_parameter.BCN_Prop;            
                end
           end
       

            clin = vertcat(clin(:));
            prop = vertcat(prop(:));   

    
            % To determine the significance for the ranking
            % Use Spearman 'rows','pairwise' to ensure that row's with NaN's in both columns are not considered in the analysis.
            [rho, pval] = corr(clin, prop,'Type','Spearman','rows','pairwise');      
          
            % Find NaN values to avoid plotting these
            idx_nan = isnan(clin) | isnan(prop);        
     
            subplot(3,1,J);
            plot(clin(~idx_nan)  , prop(~idx_nan) ,'*' )
            box off
               
            % Change fontsize
            ax = gca;
            ax.XAxis.FontSize = 14;    ax.YAxis.FontSize = 12;
               
            title(sprintf('%s, p =  %1.3f', mode{J}, pval),'FontSize',12)
            xmin = min(clin(~idx_nan));
            xmax = max(clin(~idx_nan));                
                       
            if pval < 0.05
                [P,S] = polyfit(clin(~idx_nan),prop(~idx_nan),1);
                [y_fit, ~] = polyval(P,clin(~idx_nan),S);
                
                plot(clin(~idx_nan),prop(~idx_nan),'*') 
                box off
                hold on
                
                % Plot polyfit throught data points
                plot(clin(~idx_nan),y_fit,'Color',[0.8,0.2,0.2],'LineWidth',2)
                hold on
                % Change fontsize
                ax = gca;
                ax.XAxis.FontSize = 12;    ax.YAxis.FontSize = 12;
                
                if pval < 0.001
                    title(sprintf('%s, p = <0.001, r_s = %1.3f', mode{J}, rho),'FontSize',12)
                else
                     title(sprintf('%s, p = %1.3f, r_s = %1.3f', mode{J}, pval, rho),'FontSize',12)
                end   
              
            end

            xlim([xmin xmax])
            
end

%Create one main x- and y-label
han=axes(fig,'visible','off'); 
han.YLabel.Visible='on';
han.YLabel.Position = [-0.0584    0.5000         0];
ylabel(han,{'Value SPES-prop'}, 'FontSize',15)

han.XLabel.Visible='on';
han.XLabel.Position = [ 0.5000   -0.0379         0]; 
xlabel(han,{'Value SPES-clin'}, 'FontSize',15)

% Save figure
outlabel=sprintf('All_pat_scatter.png');
path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/');
if ~exist(path, 'dir')
mkdir(path);
end
saveas(gcf,[path,outlabel],'png')    

end
         