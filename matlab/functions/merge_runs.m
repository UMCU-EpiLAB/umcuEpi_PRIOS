% This function compares the channels of the multiple runs. If these are
% equal, the runs can be merged. This function also merges these multiple
% runs of one task (for example SPES clin).

% INPUT:
% - dataBase
%   structure containing the following fields:
% - sub_label
%   string containing the name of the subject
% - ses_label
%   string containing the name of the session
% - task_label
%   string containg the name of the task
% - run_label
%   string containing the name of the run(s)
% - dataName
%   name of the eeg-file from which data is extracted
% - ccep_header
%   struct containing, among other things, the field: Fs (sample frequency)
% - ch
%   cell[channels x 1] containing the channel names 
% - tt
%   matrix[1 x timepoints] containing the time points in an epoch
% - tb_channels
%   table containing information regarding the recording channels (see BIDS
%   structure)
% - cc_stimsets
%   matrix[stim pairs x 2] containing all the electrode numbers that are
%   stimulated (C1-C2 and C2-C1 are combined)
% - cc_stimchans
%   cell[stim pairs x 2] containing all the electrode names that are
%   stimulated (C1-C2 and C2-C1 are combined)
% - cc_epoch_sorted
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined)
% - cc_epoch_sorted_avg
%   matrix[channels x stimulus pairs x samples] containing averaged
%   responses to stimulation (averaged from cc_epoch_sorted)
% - tt_epoch_sorted
%   matrix[trials x stimulus pairs x samples] containing epoched time
%   points

% OUTPUT:
% - dataBase_merge:
%   struct containing the same fields as the input struct dataBase, but now
%   the runs are merged into one. Some fields are empty, but this enables
%   automatically removal of fields before saving the data for further
%   analyses. 

function dataBase_merge = merge_runs(dataBase)

   dataBase_merge.sub_label = dataBase(1).sub_label;
   dataBase_merge.ses_label = dataBase(1).ses_label;
   dataBase_merge.task_label = dataBase(1).task_label;
   dataBase_merge.run_label = {dataBase(:).run_label};    
   dataBase_merge.dataName = dataBase(1).dataName;
   dataBase_merge.ccep_header = dataBase(1).ccep_header;
   dataBase_merge.tb_events = [];
   
   % Check whether the two tb_channels are the same, since these might be
   % different for different runs because electrodes can turn bad etc.
   for i = 1:size(dataBase,2)-1
       
        if isequal(dataBase(i).tb_channels,dataBase(i+1).tb_channels)
             dataBase_merge.tb_channels = dataBase(1).tb_channels ;           % tb_channels are equal for both SPES sessions
        
        else  % When channels are not equal, give an error. Channels should be merged manually!              
              diff_order = find(~cellfun(@isequal, dataBase(i+1).tb_channels{:,1}, dataBase(i).tb_channels{:,1}));
              diff_clin = dataBase(i).tb_channels{diff_order,1};
              diff_prop = dataBase(i+1).tb_channels{diff_order,1};
                           
              error(['For example: On row %d, %d in tb_channels, a ',...
                  'difference is found between the two runs. \nIn the 1st ',...
                  'file, it is called %s, %s, in the 2nd file %s, %s. \n',...
                  'CHECK FILES FOR OTHER DIFFERENCES!\n'], ...
                  diff_order(1), diff_order(2), diff_clin{1}, ...
                  diff_clin{2}, diff_prop{1}, diff_prop{2})

        end
   end
      
   dataBase_merge.tb_electrodes = dataBase(1).tb_electrodes;
   dataBase_merge.ch = dataBase(1).ch;
   dataBase_merge.data = [];
   dataBase_merge.Burstsup = [];
   dataBase_merge.Seizure = [];

   % Concatenate the stimpairs in the two runs
   dataBase_merge.cc_stimsets = cat(1,dataBase(:).cc_stimsets);
   dataBase_merge.cc_stimchans = cat(1,dataBase(:).cc_stimchans);         
   dataBase_merge.cc_epoch_sorted = cat(2,dataBase(:).cc_epoch_sorted);     
   dataBase_merge.cc_epoch_sorted_avg = cat(2,dataBase(:).cc_epoch_sorted_avg);         
   dataBase_merge.tt_epoch_sorted = cat(2,dataBase(:).tt_epoch_sorted);         
   dataBase_merge.tt = dataBase(1).tt;
end