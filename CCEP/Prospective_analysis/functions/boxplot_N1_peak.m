function boxplot_N1_peak(dataBase, myDataPath)
close all
clc
%% Statistics
% mode = {'Amplitude','Latency'};
new_mat = [];

% for mode_N1 = 1:size(mode,2)
    for subj = 1:size(dataBase,2)
        ccep_clin = dataBase(subj).ccep_clin;
        ccep_prop = dataBase(subj).ccep_prop;
        sub_label = dataBase(subj).sub_label;
        
        fs = 1/(size(ccep_prop.tt,2)/4);                                        % Devide tt by four because ccep_prop.tt includes 4 seconds.


%         if strcmp(mode{mode_N1},'Amplitude')
%             clin = ccep_clin.n1_peak_amplitude_check;
%             prop = ccep_prop.n1_peak_amplitude_check;

%         elseif strcmp(mode{mode_N1},'Latency')
            clin = ccep_clin.n1_peak_sample_check;
            clin = ((clin*fs)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
            prop = ccep_prop.n1_peak_sample_check;
            prop = ((prop*fs)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact

            
            
%         end
            
        % Create matrix with the clinical values in the first column and
        % the propofol values in the second
        i = 1;
        for stimp = 1:size(ccep_prop.stimpnames_avg,2)                          % For each stimpair
            for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode

            % When both clinical SPES and propofol SPES show an ER
              if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
                    new_mat(i,1) = clin(elec, stimp);            % plot the SPES-clin amp in column 1
                    new_mat(i,2) = prop(elec, stimp);          % plot the SPES-prop amp in column 2
                    i = i+1;                
              end
            end      
        end
        
        NorDisClin = lillietest(new_mat(:,1));                  % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected
        NorDisProp = lillietest(new_mat(:,2));

        if NorDisClin == 1 && NorDisProp ==1
            p = ranksum(new_mat(:,1), new_mat(:,2)) ;           % tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians
%             fprintf('The detected ERs per stimulation pair is NOT normally distributed, Wilcoxon Signed Rank test is used.\n')
            dataBase(subj).ccep_clin.p_n1 = {sprintf('%1.4f',p)};
            dataBase(subj).ccep_clin.mean_N1_lat = mean(new_mat(:,1));
            dataBase(subj).ccep_prop.p_n1 = {sprintf('%1.4f',p)};
            dataBase(subj).ccep_prop.mean_N1_lat = mean(new_mat(:,2));

        else
%             fprintf('The detected ERs per stimulation pair is normally distributed, still the Wilcoxon Signed Rank test is used.\n')
            p = ranksum(new_mat(:,1), new_mat(:,2));          % alpha default = 0.05
            dataBase(subj).ccep_clin.p_n1 = {sprintf('%1.4f',p)};
            dataBase(subj).ccep_clin.mean_N1_lat = mean(new_mat(:,1));
            dataBase(subj).ccep_prop.p_n1 = {sprintf('%1.4f',p)};
            dataBase(subj).ccep_prop.mean_N1_lat = mean(new_mat(:,2));
        end

        % Display the p value 
        if p<0.05
            fprintf('Test between the N1-Latency of %s SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is a significant difference between the two protocols \n',dataBase(subj).sub_label, p);
        else
            fprintf('Test between the N1-Latency of %s SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is NO significant difference between the two protocols \n',dataBase(subj).sub_label, p);
        end

    end
% end
    

%% Make boxPlots    
new_mat = [];  
fs = 1/(size(ccep_prop.tt,2)/4);                                        % Devide tt by four because ccep_prop.tt includes 4 seconds.
    
% for mode_N1 = 1:size(mode,2)
        
    for subj = 1:size(dataBase,2)
        ccep_clin = dataBase(subj).ccep_clin;
        ccep_prop = dataBase(subj).ccep_prop;
        sub_label = dataBase(subj).sub_label;
            
%         if strcmp(mode{mode_N1},'Amplitude')
%             clin = ccep_clin.n1_peak_amplitude_check;
%             prop = ccep_prop.n1_peak_amplitude_check;

%         elseif strcmp(mode{mode_N1},'Latency')
            clin = ccep_clin.n1_peak_sample_check;
            clin = ((clin*fs)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact
            prop = ccep_prop.n1_peak_sample_check;
            prop = ((prop*fs)-2)*1000;                                   % to convert samples to milliseconds, minus 2 becuase of the period before the stimulation artefact

%         end
        
         i = 1;
         clin_colm = 2*subj-1;                      % prealloction of the column number
         prop_colm = 2*subj;                        % prealloction of the column number

                for stimp = 1:size(ccep_prop.stimpnames_avg,2)                          % For each stimpair
                    for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode

                    % When both clinical SPES and propofol SPES show an ER
                      if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
                            new_mat(i,clin_colm) = clin(elec, stimp);            % plot the SPES-clin amp in column 1
                            new_mat(i,prop_colm) = prop(elec, stimp);          % plot the SPES-prop amp in column 2
                            i = i+1;
        %                 
                      end
                    end      
                end
                
    end
    zero_mat = find(new_mat == 0);                                          % replace zero with NaN to avoid influence on the mean
    new_mat(zero_mat) = NaN;
   

    % Create boxplot with the amplitude of SPES clin and SPES prop
    figure('Position',[205,424,1530,638])
    % columnMeans = mean(new_mat, 1, 'omitnan');
    boxplot([new_mat],'Labels',{'SPES-clin 01','SPES-prop 01','SPES-clin 02**','SPES-prop 02**','SPES-clin 03**','SPES-prop 03**','SPES-clin 04**','SPES-prop 04**','SPES-clin 05**','SPES-prop 05**','SPES-clin 06**','SPES-prop 06**'})
    title(sprintf('Latency'))
    %                 legend(findall(gca,'Tag','Box'), sprintf('mean Clin = %1.0f',columnMeans(:,1)), sprintf('mean Prop = %1.0f',columnMeans(:,2)),'Location','southoutside');

%     if strcmp(mode{mode_N1},'Amplitude')
%         ylabel('Amplitude (\muV)')
%     elseif strcmp(mode{mode_N1},'Latency')
        ylabel('Latency (milliseconds)')
%     end             

    % Save figure
    outlabel=sprintf('Latency.jpg');
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')



% end
% close Figure 2


        
%% Make scatter of the amplitude and latency
figure('Position',[302,17,938,1039])

% for mode_N1 = 1:size(mode,2)
for subj = 1:size(dataBase,2)
        
%         if strcmp(mode{mode_N1},'Amplitude')
%            axes1 = axes('Parent',figure1,'Position',[0.13,0.58,0.61,0.34]);            % For the scatterplot
           colm_clin = 2*subj-1;
           colm_prop = (2*subj);
           clin = new_mat(:,colm_clin);
           prop = new_mat(:,colm_prop);

%         elseif strcmp(mode{mode_N1},'Latency')
%            axes1 = axes('Parent',figure1,'Position',[0.13,0.11,0.63,0.34]);            % For the scatterplot 
%            clin = new_mat(:,(2*colm-1));
%            prop = new_mat(:,(2*colm));
%           end
        
%         hold(axes1,'on'); 
        subplot(size(dataBase,2),1,subj)
        scatter(clin  , prop); 
        ylabel('Propofol SPES')
        xlabel("Clinical SPES"+newline+"   ")
%%% TITLE GAAT FOUT OMDAT IK DE PVALUE RAAR OPESLAGEN HEB...%%%
%         title(fprintf('N1-%s, p =  %f, %s', mode{mode_N1}, dataBase(subj).ccep_prop.p_n1, dataBase(subj).sub_label))
        title(sprintf('N1 Latency, %s', dataBase(subj).sub_label))

        % Determine minimal and maximal value on the x-axis, used for the reference line
        % Mode amplitude and Latency are treaded differently since the
        % amplitude is often negative.

%          if strcmp(mode{mode_N1},'Amplitude')
%                 min_x = min(new_mat(:,1), [], 'all');
%                 xmin = round(min_x - 0.1*min_x, 2);
%                 max_x = max(new_mat(:,1), [], 'all');
%                 xmax = round(max_x - 0.1*max_x, 2);
%    
%                 % PLot reference line when data is significant
%                 if p < 0.05
%                     P = polyfit(new_mat(:,1), new_mat(:,2),1);
%                     X = xmin : -0.1*xmax : xmax-0.2*xmax;
%                     Y = P(1)*X + P(2);
% 
%                     hold on
%                     h=plot(X,Y);
%                     hold off
%                     h.LineWidth = 2;
%                     title(sprintf('N1-%s, %s, p = <0.05', mode{mode_N1}, dataBase(subj).sub_label))
% 
% %                     title(sprintf('%s, %s, p = <0.05', dataBase(subj).sub_label, mode{mode_N1}))
%     %                 legend(sprintf('%s',mode{mode_N1}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)      
%                     if p < 0.01
%                          title(sprintf('N1-%s, %s, p = <0.01', mode{mode_N1}, dataBase(subj).sub_label))
%                     end
% 
%                 end
%                xlim([xmin, xmax-0.2*xmax])
% 
% 
%         elseif strcmp(mode{mode_N1},'Latency')
                xmin = round(min(new_mat(:,colm_clin), [], 'all'),2);
                xmax = round(max(new_mat(:,colm_clin), [], 'all'),2);

                % PLot reference line when data is significant
                p = str2double(dataBase(subj).ccep_clin.p_n1{:});
                if p < 0.05
                    P = polyfit(new_mat(:,colm_clin), new_mat(:,colm_prop),1);
                    X = xmin : 0.1*xmax : xmax+0.2*xmax;
                    Y = P(1)*X + P(2);

                    hold on
                    h=plot(X,Y);
                    hold off
                    h.LineWidth = 2;
                    title(sprintf('N1-Latency, %s, p = <0.05', dataBase(subj).sub_label))
%                     legend(sprintf('%s',mode{mode_N1}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)      
                    if p < 0.01
                         title(sprintf('N1-Latency, %s, p = <0.01', dataBase(subj).sub_label))
                    end

                end
               xlim([xmin, xmax]) 

         end 
            


          % Save figure
        outlabel=sprintf('n1_scatter_Latency.jpg');
        path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/Scatter/');
        if ~exist(path, 'dir')
            mkdir(path);
        end
        saveas(gcf,[path,outlabel],'jpg')

end
    