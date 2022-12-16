% Observers are able to select another N1-peak than found by the detector. 
% When two observers have selected two different points, take the mean of the two observed points when the difference is <
% 0.0024 seconds. When the difference is more, than perform a visual check and
% select the correct N1-peak. 
            
function select_n1_latency(rater1, rater2, myDataPath)
      
%% load CCEPs to plot later
% this could take more than 2 minutes

files_in_folder = dir(fullfile(myDataPath.dataPath,'derivatives','CCEPs'));
[~,filename] = fileparts(rater1.dataName);
replaceRun = strfind(filename,'run');
filename = [filename(1:replaceRun-1), 'CCEP.mat'];

idx_file = contains({files_in_folder(:).name},filename);
ccep = load(fullfile(myDataPath.dataPath,'derivatives','CCEPs',files_in_folder(idx_file).name)); 

%% correct N1-latencies

% preallocation
n1_peak_sample_check_comb = NaN(size(rater1.n1_peak_amplitude_check));

for stimp = 1:size(rater1.n1_peak_sample_check,2)
    for elec = 1:size(rater1.n1_peak_sample_check,1)
        if ~isnan(rater1.n1_peak_sample_check(elec,stimp)) && ~isnan(rater2.n1_peak_sample_check(elec,stimp))
            
            if abs(rater1.n1_peak_sample_check(elec,stimp) - rater2.n1_peak_sample_check(elec,stimp)) > 5 
                
                % plot figure to determine whether observer 1 or 2 was right
                H = figure(1);
                H.Units = 'normalized';
                H.Position = [0.14,0.0625,0.77,0.7];
    
                plot(ccep.tt, squeeze(ccep.cc_epoch_sorted_reref_avg(elec,stimp,:)),'k','linewidth',2);  % plot the rereference signal in a solid line
                hold on
                
                plot(ccep.tt(rater1.n1_peak_sample_check(elec,stimp)), rater1.n1_peak_amplitude_check(elec,stimp),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
                plot(ccep.tt(rater2.n1_peak_sample_check(elec,stimp)), rater2.n1_peak_amplitude_check(elec,stimp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',4)
                plot(ccep.tt,squeeze(ccep.cc_epoch_sorted_reref(elec,stimp,:,:)),'r:');                % plot the 10 separate stimulations

                % Create a patch for the 0-9 ms interval post-stimulation
                % in which no physiological activity can be measured.
                patch([0 19/2048 19/2048 0],[-800 -800 750 750],[0.6,0.2,0.2],'EdgeAlpha',0)
                alpha(0.2)
                
                hold off
    
                xlim([-0.05 0.2])
                ymin = floor(1.1*rater1.n1_peak_amplitude(elec,stimp));
                ylim([ymin 600])
                xlabel('Time(s)')
                ylabel('Potential (\muV)')
                title(sprintf('Electrode %s, stimulating %s-%s',ccep.ch{elec},ccep.cc_stimchans{stimp,:}))
                legend(' ','Rater1','Rater2',' ')

                % Decide which of the two N1-peaks is the correct one
                answer = input('Which N1-peak is correct selected? Rater1(red) or Rater2(blue) [1/2]: ');

                if answer == 1
                    n1_peak_sample_check_comb(elec,stimp) = rater1.n1_peak_sample_check(elec,stimp);

                elseif answer == 2 
                    n1_peak_sample_check_comb(elec,stimp) = rater2.n1_peak_sample_check(elec,stimp);

                end

            else % when latencies for both observers were (almost) equal  
                n1_peak_sample_check_comb(elec,stimp) = round(mean([rater1.n1_peak_sample_check(elec,stimp),rater2.n1_peak_sample_check(elec,stimp)]));

            end
        end
    end
end

%% save file
% save struct in which visual checks of the N1-latencies of both observers are combined

raterComb.n1_peak_sample_check  = n1_peak_sample_check_comb;
raterComb.ch                    = rater1.ch;
raterComb.dataName              = rater1.dataName;
raterComb.n1_peak_sample        = rater1.n1_peak_sample;
raterComb.cc_stimchans          = rater1.cc_stimchans;
raterComb.cc_stimsets           = rater1.cc_stimsets;
raterComb.tt                    = rater1.tt;

filefolder = fullfile(myDataPath.CCEPpath,'checkedN1s');
[~,filename] = fileparts(rater1.dataName);
replaceRun = strfind(filename,'run');
filename = [filename(1:replaceRun-1), 'N1sChecked_comb.mat'];

if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

save(fullfile(filefolder,filename),'-struct','raterComb','-v7.3');

fprintf('N1-latencies are checked and saved as: \n%s \n',fullfile(filefolder,filename))

end % end function



