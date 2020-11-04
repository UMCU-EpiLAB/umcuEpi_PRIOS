function  scatter_networkPar(dataBase, myDataPath)

mode = {'ERs per stimulation pair','Indegree','Outdegree','BC'};

for J = 1:size(mode,2)
    
    figure('Position',[302,17,1224,1039])
    
    for i = 1:size(dataBase,2)
        
        if strcmp(mode{J},'ERs per stimulation pair')
            par10 = dataBase(i).agreement_parameter.ERs_stimp10;
            par2 = dataBase(i).agreement_parameter.ERs_stimp2;
            pval = dataBase(i).statistics.p_stimp;
            rho = dataBase(i).statistics.rho_stimp;
            
        elseif strcmp(mode{J},'Indegree')
            par10 = dataBase(i).agreement_parameter.indegreeN_10;
            par2 = dataBase(i).agreement_parameter.indegreeN_2;
            pval = dataBase(i).statistics.p_indegree;
            rho = dataBase(i).statistics.rho_indegree;
            
        elseif strcmp(mode{J},'Outdegree')
            par10 = dataBase(i).agreement_parameter.outdegreeN_10;
            par2 = dataBase(i).agreement_parameter.outdegreeN_2;
            pval = dataBase(i).statistics.p_outdegree;
            rho = dataBase(i).statistics.rho_outdegree;
            
        elseif strcmp(mode{J},'BC')
            par10 = dataBase(i).agreement_parameter.BCN_10;
            par2 = dataBase(i).agreement_parameter.BCN_2;
            pval = dataBase(i).statistics.p_BC;
            rho = dataBase(i).statistics.rho_BC;
        end
        
        subplot(size(dataBase,2),1,i)
        scatter(par10  , par2  )
        ylabel('2 stimuli setting')
        xlabel("10 stimuli setting"+newline+"   ")
        str_main = sprintf('%s', mode{J});
        sgtitle(str_main)
        title(sprintf('%s, p =  %1.3f', dataBase(i).sub_label, pval))
        legend(sprintf('%s',mode{J}),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)
        xmin = 0;
        xmax = round(max(par10)+0.1*max(par10),2);
        
        if pval < 0.05
%             num_at_zero  = numel(par2(par10 ==0))+1;
%             intercept_zero = sum(par2(par10 ==0))/num_at_zero;           % Find the mean value for the interception point with the zero line for the reference line
%             h = refline(rho,intercept_zero);
            idx_nan = isnan(par10) | isnan(par2);
            P = polyfit(par10(~idx_nan),par2(~idx_nan),1);
            X = xmin:0.1*xmax:xmax+0.2*xmax;
            Y = P(1)*X + P(2);
            
            %             intercept_zero = mean(par2(par10 ==0));           % Find the mean value for the interception point with the zero line for the reference line
            %             h = refline(rho,intercept_zero);
            hold on
            h=plot(X,Y);
            hold off
            h.LineWidth = 2;
            title(sprintf('%s, p = <0.05', dataBase(i).sub_label))
            legend(sprintf('%s',mode{J}), sprintf('r_s = %1.3f',rho  ),'Location','EastOutside','Orientation','vertical','Box','off','FontSize',12)
            if pval < 0.01
                title(sprintf('%s, p = <0.01', dataBase(i).sub_label))
            end
            
        end
        xlim([xmin xmax])
        
    end
    
    % Save figure
    outlabel=sprintf('All_pat_scatter_%s.jpg',mode{J});
    path = fullfile(myDataPath.CCEPpath,'Visualise_agreement/Scatter/');
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'jpg')
    
end

end


