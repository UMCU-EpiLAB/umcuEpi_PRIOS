function plot_all_ccep(ccep_clin, ccep_prop, myDataPath)

tt = ccep_clin.tt;    
set(groot,'defaultFigureVisible','on') % 'on' to turn figures showing on, 'off' to not show the figures.  

% When the status of the channel is bad, not visualise the response
bad_sig_clin(:,1) = find(contains(ccep_clin.tb_channels.status, {'bad'}));  
ccep_clin.cc_epoch_sorted_avg(bad_sig_clin,:,:) = NaN;

bad_sig_prop(:,1) = find(contains(ccep_prop.tb_channels.status, {'bad'}));
ccep_prop.cc_epoch_sorted_avg(bad_sig_prop,:,:) = NaN;

%%% WERKT NIET ALS DE EEN EEN MINDER STIMPAIREN HEEFT DAN DE ANDER (30-9)
% Find rows (stimpair) which are not equal for both runs
if size(ccep_clin.stimpnames_avg,2) ~= size(ccep_prop.stimpnames_avg,2)
   diff_row = find(~ismember(ccep_clin.stimpnames_avg, ccep_prop.stimpnames_avg));      % Find stimpairs which are in clin-SPES but not in prop-SPES
   warning('Clinical SPES contains more stimpairs than Prop-SPES')
   
   ccep_clin.stimpnames_avg(diff_row) = [];
end

% Plot all averaged responses to the stimuli of one stimulation pair
% Left the Clinical SPES, right the Propofol SPES
% Only plot the stimulation pairs which are stimulated in both protocols
for stimp = 1:length(ccep_clin.stimpnames_avg)            % For each stimulation pair                              
    Stimpnm_clin = ccep_clin.stimpnames_avg{stimp};     
    Stimpnm_prop = ccep_prop.stimpnames_avg{stimp};     

    % When the stimulation pairs are equal for clin and prop, plot them
    % in one figure to compare
    if Stimpnm_clin == Stimpnm_prop
        % Clinical SPES
        % Electrodes in the stimulation pair
        stim_elec = ccep_clin.cc_stimsets_avg(stimp,:);
       
        % Remove electrodes which are stimulated
        plot_ccep = ccep_clin.cc_epoch_sorted_avg(:,:,:);
        plot_ccep(stim_elec,:,:) = NaN;

        plot_ccep_clin = squeeze(plot_ccep(:,stimp,:));
        plot_ccep_clin(:,tt>-0.01 & tt<0.02) = NaN;                 % remove stimulation artefact for easier visual check

        % Propofol SPES
        % Electrodes in the stimulation pair
        stim_elec_prop = ccep_prop.cc_stimsets_avg(stimp,:);
        
        % Remove electrodes which are stimulated
        plot_ccep_prop = ccep_prop.cc_epoch_sorted_avg(:,:,:);
        plot_ccep_prop(stim_elec_prop,:,:) = NaN;

        plot_ccep_prop = squeeze(plot_ccep_prop(:,stimp,:));
        plot_ccep_prop(:,tt>-0.01 & tt<0.02) = NaN;
        
        
        % Plot Clinical SPES (left) and Propofol SPES (right)
        figure('Position',[673,21,1227,1041])
        sub1 =  subplot(1,2,1);
        sub1.Position = [0.07,0.11,0.41,0.82];
        plot(tt, plot_ccep_clin' + [0:1000:size(plot_ccep_clin,1)*1000-1])

        str = sprintf('%s Clinical SPES', Stimpnm_clin);
        title(str)
        set(gca,'YTick',[0:1000:size(plot_ccep_clin,1)*1000-1],'YTickLabel',ccep_clin.ch) ;                  
        xlim([-.2 1.5])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 


        % plot propofol SPES
        sub2 = subplot(1,2,2);
        sub2.Position = [0.57,0.11,0.41,0.82];
        plot(tt, plot_ccep_prop' + [0:1000:size(plot_ccep_prop,1)*1000-1])

        str = sprintf('%s Propofol SPES', Stimpnm_prop);
        title(str)
        set(gca,'YTick',[0:1000:size(plot_ccep_prop,1)*1000-1],'YTickLabel',ccep_prop.ch) ;                  
        xlim([-.2 1.5])
        ylabel('All stimuli of this stimulation pair' )
        xlabel('time (s)') 



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
    
end
close all
end
            

          
