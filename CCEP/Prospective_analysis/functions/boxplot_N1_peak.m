function boxplot_N1_peak(dataBase, myDataPath)
close all
    
figure1 = figure('Position',[349,383,891,679]);         % This is for the scatter plot 
new_mat = [];   
    

    mode = {'Amplitude','Latency'};

for mode_N1 = 1:size(mode,2)
        
    for subj = 1:size(dataBase,2)
        ccep_clin = dataBase(subj).ccep_clin;
        ccep_prop = dataBase(subj).ccep_prop;
        sub_label = dataBase(subj).sub_label;
            
        if strcmp(mode{mode_N1},'Amplitude')
            clin = ccep_clin.n1_peak_amplitude_check;
            prop = ccep_prop.n1_peak_amplitude_check;
            axes1 = axes('Parent',figure1,'Position',[0.13,0.58,0.61,0.34]);            % For the scatterplot


        elseif strcmp(mode{mode_N1},'Latency')
            clin = ccep_clin.n1_peak_sample_check;
            prop = ccep_prop.n1_peak_sample_check;
            axes1 = axes('Parent',figure1,'Position',[0.13,0.11,0.63,0.34]);            % For the scatterplot      

        end

%         % Find the number of times a ER is found in SPES-clin and SPES prop
%         ER_clin = find(~isnan(clin)==1 );
%         ER_prop = find(~isnan(prop)==1);
%         ERs_both = find(ismember(ER_clin, ER_prop,'rows'));
        

        i = 1;
%%% DIT WERKT NIET VANAF 4

        if subj > 1
            subj = subj+1;
        end
        
        if size(ccep_prop.stimpnames_avg,2) == size(ccep_clin.stimpnames_avg,2)         % Both runs must have the same number of stimpairs
            if size(ccep_prop.ch,1) == size(ccep_clin.ch,1)                             % Both runs must have the same number of electrodes
                for stimp = 1:size(ccep_prop.stimpnames_avg,2)                          % For each stimpair
                    for elec = 1:size(ccep_prop.ch,1)                                   % For each electrode

                    % When both clinical SPES and propofol SPES show an ER
                      if ~isnan(clin(elec, stimp)) &&  ~isnan(prop(elec, stimp)) 
                            new_mat(i,subj) = clin(elec, stimp);          % plot the SPES-clin amp in column 1
                            new_mat(i,subj+1) = prop(elec, stimp);          % plot the SPES-prop amp in column 2
                            i = i+1;
        %                 
                      end
                    end      
                end
            end 
        end
    end
    zero_mat = find(new_mat == 0);
    new_mat(zero_mat) = NaN;
    

end
%                 mean_clin = mean(new_mat(:,1));
%                 mean_prop = mean(new_mat(:,2));
                % Create boxplot with the amplitude of SPES clin and SPES prop
                figure('Position',[348,362,892,700])
                boxplot([new_mat],'Labels',{'SPES-clin 01','SPES-prop 01','SPES-clin 02','SPES-prop 02','SPES-clin 03','SPES-prop 03','SPES-clin 04','SPES-prop 04','SPES-clin 05','SPES-prop 05','SPES-clin 06','SPES-prop 06'})
                title(sprintf('%s, %s',sub_label{1},mode{mode_N1}))
%                 legend(findall(gca,'Tag','Box'), sprintf('mean Clin = %1.0f',mean_clin), sprintf('mean Prop = %1.0f',mean_prop),'Location','southeastoutside');

                if strcmp(mode{mode_N1},'Amplitude')
                    ylabel('Amplitude (\muV)')
                elseif strcmp(mode{mode_N1},'Latency')
                    ylabel('Latency (samples)')
                end
%             else
%                 warning('CCEP Clin and CCEP Prop do not have the same number of electrodes')
%             end
%         else
%             warning('CCEP Clin and CCEP Prop do not have the same number of stimulation pairs')
%         end

        % Save figure
        outlabel=sprintf('sub-%s_%s.jpg',...
            sub_label{1},mode{mode_N1});
        path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/');
        if ~exist(path, 'dir')
            mkdir(path);
        end
        saveas(gcf,[path,outlabel],'jpg')

        close Figure 2


    %% Make scatter of the amplitude and latency

        NorDisClin = lillietest(new_mat(:,1));                  % null hypothesis that x is normally distributed, results in 1 when the null hypothesis is rejected 
        NorDisProp = lillietest(new_mat(:,2));

        if NorDisClin == 1 && NorDisProp ==1
            p = ranksum(new_mat(:,1), new_mat(:,2)) ;           % tests the null hypothesis that data in x and y are samples from continuous distributions with equal medians
            fprintf('The detected ERs per stimulation pair is NOT normally distributed, Wilcoxon Signed Rank test is used.\n')

        else
            fprintf('The detected ERs per stimulation pair is normally distributed, still the Wilcoxon Signed Rank test is used.\n')
            p = ranksum(new_mat(:,1), new_mat(:,2));          % alpha default = 0.05

        end

        % Display the p value 
        if p<0.05
            fprintf('Test between the N1-%s of SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is a significant difference between the two protocols for %s \n', mode{mode_N1},p, sub_label{1});
        else
            fprintf('Test between the N1-%s of SPES-clin and SPES-prop gives p-value = %1.4f. This means that there is NO significant difference between the two protocols for %s \n',mode{mode_N1}, p, sub_label{1});
        end


        hold(axes1,'on'); 
        scatter(new_mat(:,1)  , new_mat(:,2), 'ok','Parent',axes1); 
        ylabel('Propofol SPES')
        xlabel("Clinical SPES"+newline+"   ")
        title(sprintf('N1-%s, p =  %1.3f, %s', mode{mode_N1}, p, sub_label{1}))
    %     legend(sprintf('%s',mode{mode_N1}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)      


        % Determine minimal and maximal value on the x-axis, used for the reference line
        % Mode amplitude and Latency are treaded differently since the
        % amplitude is often negative.

         if strcmp(mode{mode_N1},'Amplitude')
                min_x = min(new_mat(:,1), [], 'all');
                xmin = round(min_x - 0.1*min_x, 2);
                max_x = max(new_mat(:,1), [], 'all');
                xmax = round(max_x - 0.1*max_x, 2);
    %             set(gca,'xtick',[])
    %             set(gca,'ytick',[])

                % PLot reference line when data is significant
                if p < 0.05
                    P = polyfit(new_mat(:,1), new_mat(:,2),1);
                    X = xmin : -0.1*xmax : xmax-0.2*xmax;
                    Y = P(1)*X + P(2);

                    hold on
                    h=plot(X,Y);
                    hold off
                    h.LineWidth = 2;
                    title(sprintf('%s, %s, p = <0.05', sub_label{1}, mode{mode_N1}))
    %                 legend(sprintf('%s',mode{mode_N1}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)      
                    if p < 0.01
                        title(sprintf('%s, %s, p = <0.01', sub_label{1},mode{mode_N1}))
                    end

                end
               xlim([xmin, xmax-0.2*xmax])


        elseif strcmp(mode{mode_N1},'Latency')
                xmin = round(min(new_mat(:,1), [], 'all'),2);

                xmax = round(max(new_mat(:,1), [], 'all'),2);
    %             set(gca,'xtick',[])
    %             set(gca,'ytick',[])

                % PLot reference line when data is significant
                if p < 0.05
                    P = polyfit(new_mat(:,1), new_mat(:,2),1);
                    X = xmin : 0.1*xmax : xmax+0.2*xmax;
                    Y = P(1)*X + P(2);

                    hold on
                    h=plot(X,Y);
                    hold off
                    h.LineWidth = 2;
                    title(sprintf('%s, %s, p = <0.05', sub_label{1}, mode{mode_N1}))
    %                 legend(sprintf('%s',mode{mode_N1}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)      
                    if p < 0.01
                        title(sprintf('%s, %s, p = <0.01', sub_label{1}, mode{mode_N1}))
                    end

                end
               xlim([xmin, xmax]) 

         end 

%     end 
          % Save figure
        outlabel=sprintf('sub-%s_n1_peak_scatter.jpg',sub_label{1});
        path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/N1_compare/Scatter/');
        if ~exist(path, 'dir')
            mkdir(path);
        end
        saveas(gcf,[path,outlabel],'jpg')
% end
    end