function plot_all_ccep(dataBase, myDataPath, LocOnes, stimchans, dif_mat)

tic;
ccep_plot = zeros(10,8192);
tt = dataBase.tt;    
LocaOnes = LocOnes{:,:}; 

% ER in 10 stims = 1, ER in 2 stims = -1 (dif_mat = 10stim - 2stim --> ER in 10stim (1) - non-ER in 2 stim(0) = 1)
[ER_in10(:,2), ER_in10(:,1)] = find(dif_mat == 1);            % [stimpairs electrodes]

for i = 1:size(ER_in10(:,2))                            % For the number of ones detected
    ER_in10st(i,1) = stimchans(ER_in10(i,1))';
    ER_in10st(i,2) = dataBase.ch(ER_in10(i,2));    
end


set(groot,'defaultFigureVisible','on') % 'on' to turn figures showing on, 'off' to not show the figures.  

[indivstimp,~,stimprow] = unique(sort(dataBase.cc_stimsets,2),'rows');


for stimp = 1:length(indivstimp)                            %1:max(indivstimp)
           
        Stimpnm = stimchans{stimp};     
        stimpnm = find(stimprow == stimp);

       for elec = 1:size(dataBase.cc_epoch_sorted_avg,1)                     % for every electrode
            elecnm = dataBase.ch{elec};
            
            for i = 1:size(LocaOnes)     
                if ismember({Stimpnm}, LocaOnes{i,1}) && ismember(dataBase.ch{elec}, LocaOnes(i,2)) 

                    test = squeeze(dataBase.cc_epoch_sorted(elec,:,stimpnm,:));
                    test2 = reshape(test, size(test,1)*size(test,2), size(test,3));
                    test2(:,tt>-0.01 & tt<0.02) = NaN;
                
                    figure('Position',[1100 0 700 700])
                    plot(tt,test2'+[0:500:size(test2,1)*500-1]);                    
                    str = sprintf('ER in 2 stims. Stimulation pair %s for electrode %s', Stimpnm, elecnm);
                    title(str)
                    hold on
                    
                    for j = 1:length(ER_in10st)
                        if ismember({Stimpnm}, ER_in10st{j,1}) && ismember(dataBase.ch{elec}, ER_in10st(j,2)) 
                            plot(tt, test2'+[0:500:size(test2,1)*500-1]);           
                            str = sprintf('ER in 10 stims. Stimulation pair %s for electrode %s', Stimpnm, elecnm);
                            title(str)
                        end
                    end
                    

                    set(gca,'YTick',500*(0:size(ccep_plot)-1))
                    
                    xlim([-.2 1.5])
                    ylabel('All stimuli of this stimulation pair' )
                    xlabel('time (s)') 
                    hold off

        % Save the figures
            if dataBase.save_fig==1
                % create folder to save figures
                if ~ exist(fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm),'dir')

                    mkdir(fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm));
                end

                % filename
                figureName = fullfile(myDataPath.CCEPpath,'ccep_figures',dataBase.sub_label,Stimpnm,...
                    [dataBase.sub_label '_stimp_' Stimpnm '_elec_' elecnm ]);
                set(gcf,'PaperPositionMode','auto');
                print('-dpng','-r300',figureName);
            else
                pause
            end
            
        end
     end
            
            
       end
        close all          
   end
    toc
end

          
