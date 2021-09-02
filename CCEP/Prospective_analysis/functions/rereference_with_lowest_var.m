
function dataBase = rereference_with_lowest_var(dataBase,cfg)
% The signals can be re-reference by a reference.
% The reference calculatede below is the median of the signals that have a
% post- and pre-stim variance that is smaller than the pre-stim variance of
% the median signal of all signals (CAR). At least 5% of the signals have
% to be averaged in the reference signal. 

% Lowest variance is used to avoid introducing CCEP's in signals that
% originally do not have a CCEP.

period_preStim = 3072:3891;                             % [-0.5: -0.1] pre-stim period, fs = 2048 --> [4096-1024 : 4096-205]
period_postStim = 4116:4301;                            % [0.01: 0.1] post-stim period --> [4096+20 : 4096+205]

% for subj = 1:size(dataBase,2)
%     for run = 1:size(dataBase.run_label,2)          % not required because of merge_runs.m
       

        if strcmp(cfg.reref,'y')
            for stimp = 1:size(dataBase.cc_epoch_sorted_select,2)                  % For each stimulation pair. 
                                
                % Exclude bad channels and stimulated channels by making it NaN's
                data_all = dataBase.cc_epoch_sorted_select;
                bad_and_stimp = [find(strcmp(dataBase.tb_channels.status_description,'noisy (visual assessment)')); dataBase.cc_stimsets_avg(stimp,:)'];
                data_all(bad_and_stimp,:,:,:) = NaN;
                        
                for numstim = 1:size(data_all,3)                                                % for each of the trials of one stimulus pair

                    % Determine median of all signals except BAD and stimulated channels (CAR)
                    raw_data_stim = squeeze(data_all(:,stimp,numstim,:));                  
                    CAR_raw_data = median(raw_data_stim,'omitnan'); 


                    % Determine variance of CAR signal PRE-stim
                    var_preStimCAR = var(CAR_raw_data(:, period_preStim),'omitnan');


                    % Determine for each separate signal if variance
                    % POST-stim is lower than variance of CAR signal
                    % var(x,0 or 1) for normalization (Normalization is a
                    % good technique to use when you do not know the distribution of your data)
                    % var(x,1,1) = var for columns, var(x,1,2) = var for rows
                    var_postStim_signal = var(squeeze(data_all(:,stimp,numstim, period_postStim)),1,2,'omitnan');


                    % Determine for each separate signal if variance
                    % pre-stim is lower than variance of CAR signal
                    var_preStim_signal = var(squeeze(data_all(:,stimp,numstim, period_preStim )),1,2,'omitnan');


                    % Keep the signals that have a variance post- AND pre-stim that is lower than the variance of the CAR signal
                    % Determine the median --> this is your reference signal 
                    post_lower = find(var_postStim_signal < var_preStimCAR);
                    pre_lower = find(var_preStim_signal < var_preStimCAR);
                    ref_keep = post_lower(ismember(post_lower,pre_lower));

                    ref = median(squeeze(data_all(ref_keep, stimp,numstim, :)));

                    % If there are less than 10% of the signals present in the reference.
                    % Then add signals with lowest variance POST-stim.
                    if size(ref_keep,1) < ceil(0.1*size(data_all,1))

                        [~,idx_var_post] = sort(var_postStim_signal,'ascend');                              % NAN data receives the highest scoring
                        nmb_extra_signals = round(0.1*size(raw_data_stim,1))- size(ref_keep,1);            % calculate the number of extra signals to obtain 5% of the signals in the reference

                        % If the lowest 10% variance is already in the ref_keep, then take the number of extra required
                        % signals + the number signals that were already present in the ref_keep
                        if any(ismember(ref_keep, idx_var_post(1:(size(ref_keep,1)+nmb_extra_signals))))
                           nmb_extra_signals = nmb_extra_signals + sum(ismember(ref_keep, idx_var_post(1:(size(ref_keep,1)+nmb_extra_signals))));   
                        end

                        ref_keep2 = unique([ref_keep; idx_var_post(1:nmb_extra_signals)]);

                        ref = median(squeeze(data_all(ref_keep2, stimp,numstim, :)));
                    end

                    dataBase.ref(:,stimp,numstim,:) = ref;


                    % Re-reference the individual trials. Make sure to take the original signal to make
                    % sure to also rereference the bad channels and stimulated channels
                    dataBase.cc_epoch_sorted_select_reref(:,stimp,numstim,:) = squeeze(dataBase.cc_epoch_sorted_select(:,stimp,numstim,:)) - ref;
                end

            end
            

        end
       
%     end
% end

% Take mean of all numstims of the re-referenced signals per stimulation pair
dataBase.cc_epoch_sorted_select_reref_avg = squeeze(mean(dataBase.cc_epoch_sorted_select_reref,3,'omitnan')); 

fprintf('%s has been re-referenced. \n',dataBase.sub_label)


end

%% Old script to visualise result in figure

% Espcially visualise the stims with artefact periods.
% 
% 
% tt = (1:cfg.epoch_length  *2048)/2048 - cfg.epoch_prestim;
% 
% chan = 42;
% stimp = 35;
% numstim = 1;
% 
% figure()
% plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,numstim,:))'-300,'color',[1,0,0,0.3])
% hold on
% plot(tt,squeeze(dataBase.cc_epoch_sorted_select_reref(chan,stimp,numstim,:))'-300,'color',[0.58,0.5,0.99])
% 
% plot(tt,squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,numstim,:))'+200,'color',[1,0,0,0.3])
% plot(tt, (squeeze(dataBase.cc_epoch_sorted_select(chan,stimp,numstim,:))'- median(squeeze(data_all(:,stimp,numstim,:)),'omitnan')) +200,'color','m')
% 
% plot(tt, squeeze(dataBase.ref(:,stimp,numstim,:))' + 600)
% plot(tt, median(squeeze(data_all(:,stimp,numstim,:)),'omitnan') + 800)
% xlim([-0.5, 1.5])
% ylim([-800 1100])
% title(sprintf('Electrode %s, stimulating %s-%s',...
% dataBase.ch{chan},...
% dataBase.cc_stimchans_avg{stimp,1},...
% dataBase.cc_stimchans_avg{stimp,2}))
% xlabel('time (sec)')
% 
% 
% % Create a patch for the -1/5:9 ms interval in which no
% % physiological activity can be measured.
% patch([0 0.009 0.009 0],[-1000 -1000 size(gcf,1)*1000-1 size(gcf,1)*1000-1],[0.6,0.2,0.2],'EdgeAlpha',0)
% alpha(0.2)
% 
% legend('original','rereferenced with ref','original','rereferenced with CAR','ref','CAR')
