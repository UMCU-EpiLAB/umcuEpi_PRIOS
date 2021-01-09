function plot_all_ccep(ccep_clin, ccep_prop, myDataPath)
% Plot the average signal of all electrodes per stimulation pair
% Create figures with a plot per stimulation pair with the averaged response per electrode 
% Easy compare between clinical-SPES and propofol-SPES

tt = ccep_clin.tt;    
set(groot,'defaultFigureVisible','on')                                                  % 'on' to turn figures showing on, 'off' to not show the figures.  

% When the status of the channel is bad, not visualise the response
bad_sig_clin(:,1) = find(contains(ccep_clin.tb_channels.status, {'bad'}));  
ccep_clin.cc_epoch_sorted_avg(bad_sig_clin,:,:) = NaN;

bad_sig_prop(:,1) = find(contains(ccep_prop.tb_channels.status, {'bad'}));
ccep_prop.cc_epoch_sorted_avg(bad_sig_prop,:,:) = NaN;


% Find rows (stimpair) which are not equal for both runs  
if sum(~ismember(ccep_clin.stimpnames_avg, ccep_prop.stimpnames_avg)) ~= 0
    
    diff_row = find(~ismember(ccep_clin.stimpnames_avg, ccep_prop.stimpnames_avg));     % Find stimpairs which are in clin-SPES but not in prop-SPES
    names = ccep_clin.stimpnames_avg(diff_row);                                         % Find names which are different for both protocols (could be a typo)
    stringsz = [repmat('%s, ',1,size(names,2)-1),'%s '];                                % Create array in the lenght of the number of diff stimpairs
    sprintf(['Runs present in SPESclin and NOT in SPESprop: \n' stringsz '\n'],names{:})  % stimpairs not in both runs
      
    % The above selected incorrect stimulation pair is removed from the
    % matrices
    ccep_clin.stimpnames_avg(diff_row) = [];
    ccep_clin.cc_epoch_sorted_avg(:,diff_row,:) = [];
    ccep_clin.cc_stimsets_avg(diff_row,:) = [];  
    
    % Find stimpairs which are in SPES-prop but not in SPES-clin
    diff_row_prop = find(~ismember(ccep_prop.stimpnames_avg, ccep_clin.stimpnames_avg));
    names_prop = ccep_prop.stimpnames_avg(diff_row_prop);
    stringsz = [repmat('%s, ',1,size(names_prop,2)-1),'%s '];
    sprintf(['Runs present in SPESprop and NOT in SPESclin: \n' stringsz '\n'],names_prop{:})  % stimpairs not in both runs
    ccep_prop.stimpnames_avg(diff_row_prop) = [];
    ccep_prop.cc_epoch_sorted_avg(:,diff_row_prop,:) = [];
    ccep_prop.cc_stimsets_avg(diff_row_prop,:) = [];

end


% Plot all averaged responses to the stimuli of one stimulation pair
% Left the Clinical SPES, right the Propofol SPES
% Only plot the stimulation pairs which are stimulated in both protocols

%%% GAAT NOG FOUT VOOR PRIOS05, STIMPAIR WORDT AFGEBEELD DUS GAAT IETS NIET
%%% GOED IN HET STUK, UITZOEKEN WELK STIMPAAR HET IS.

for stimp = 1:length(ccep_clin.stimpnames_avg)            % For each stimulation pair should now be similar for prop and clin, see above                            
    Stimpnm_clin = ccep_clin.stimpnames_avg{stimp};     
    Stimpnm_prop = ccep_prop.stimpnames_avg{stimp};     

        %%% Clinical SPES
        % Electrodes in the stimulation pair
        stim_elec = ccep_clin.cc_stimsets_avg(stimp,:);
       
        % Remove electrodes which are stimulated
        plot_ccep = ccep_clin.cc_epoch_sorted_avg(:,:,:);
        plot_ccep(stim_elec,:,:) = NaN;
        plot_ccep_clin = squeeze(plot_ccep(:,stimp,:));
            
        %%% Propofol SPES
        % Electrodes in the stimulation pair
        stim_elec_prop = ccep_prop.cc_stimsets_avg(stimp,:);
        
        % Remove electrodes which are stimulated
        plot_ccep_prop = ccep_prop.cc_epoch_sorted_avg(:,:,:);
        plot_ccep_prop(stim_elec_prop,:,:) = NaN;
        plot_ccep_prop = squeeze(plot_ccep_prop(:,stimp,:));
        
        % Plot Clinical SPES (left) and Propofol SPES (right)
        figure('Position',[673,21,1227,1041])
        sub1 =  subplot(1,2,1);
        sub1.Position = [0.07,0.11,0.41,0.82];
        plot(tt, plot_ccep_clin' + (0:1000:size(plot_ccep_clin,1)*1000-1))

        str = sprintf('%s Clinical SPES', Stimpnm_clin);
        title(str)
        set(gca,'YTick',(0:1000:size(plot_ccep_clin,1)*1000-1),'YTickLabel',ccep_clin.ch) ;                  
        xlim([-.2 1.5])
        ylim([-1000 size(plot_ccep_clin,1)*1000-1])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 
        
        % Create a patch for the -1/5:9 ms interval in which no
        % physiological activity can be measured.
        patch([0 0.009 0.009 0],[-1000 -1000 size(plot_ccep_clin,1)*1000-1 size(plot_ccep_clin,1)*1000-1],[0.6,0.2,0.2],'EdgeAlpha',0)
        alpha(0.2)        

        % plot propofol SPES
        sub2 = subplot(1,2,2);
        sub2.Position = [0.57,0.11,0.41,0.82];
        plot(tt, plot_ccep_prop' + (0:1000:size(plot_ccep_prop,1)*1000-1))

        str = sprintf('%s Propofol SPES', Stimpnm_prop);
        title(str)
        set(gca,'YTick',(0:1000:size(plot_ccep_prop,1)*1000-1),'YTickLabel',ccep_prop.ch) ;                  
        xlim([-.2 1.5])
        ylim([-1000 size(plot_ccep_clin,1)*1000-1])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 
        
        % Create a patch for the -1/5:9 ms interval in which no
        % physiological activity can be measured.
        patch([0 0.009 0.009 0],[-1000 -1000 size(plot_ccep_clin,1)*1000-1 size(plot_ccep_clin,1)*1000-1],[0.6,0.2,0.2],'EdgeAlpha',0)
        alpha(0.2)
        
        
       % Save the figures
            if ccep_clin.save_fig=='y'
                % create folder to save figures
                if ~ exist(fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_clin.sub_label),'dir')

                    mkdir(fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_clin.sub_label));
                end

                % filename
                figureName = fullfile(myDataPath.CCEPpath,'all_CCEP_StimP',ccep_clin.sub_label,...
                    [ccep_clin.sub_label '_stimp_' Stimpnm_clin ]);
                set(gcf,'PaperPositionMode','auto');
                print('-dpng','-r300',figureName);
            else
                pause
            end
        
end
close all
end
            

          
