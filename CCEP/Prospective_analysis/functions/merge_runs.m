function dataBase_merge = merge_runs(dataBase)

   dataBase_merge.sub_label = dataBase(1).sub_label;
   dataBase_merge.ses_label = dataBase(1).ses_label;
   dataBase_merge.task_label = dataBase(1).task_label;
   dataBase_merge.run_label = {dataBase(:).run_label};    
   dataBase_merge.dataName = dataBase(1).dataName;
   dataBase_merge.ccep_header = dataBase(1).ccep_header;
   dataBase_merge.ch = dataBase(1).ch;
   dataBase_merge.max_stim = dataBase(1).max_stim;
   dataBase_merge.tt = dataBase(1).tt;
   
   % Check whether the two tb_channels are the same, since these might be
   % different for different runs because electrodes can turn bad etc.
   for i = 1:size(dataBase,1)-1
       
        if isequal(dataBase(i).tb_channels,dataBase(i+1).tb_channels)
             dataBase_merge.tb_channels = dataBase(1).tb_channels ;           % tb_channels are equal for both SPES sessions
        else
              diff_chan = setdiff(dataBase(i).tb_channels, dataBase(i+1).tb_channels);
              warning('These channels are different in tb_channels CHECK the diference %s \n', diff_chan.name{:})
        
             % The difference is probably little, therefore the tb_channels in
             % the merge dataBase is set anyway to make sure it is not forgotten
             dataBase_merge.tb_channels = dataBase(1).tb_channels ;
        end
   end
      
   % Concatenate the stimpairs in the two runs
   dataBase_merge.cc_stimsets_all = cat(1,dataBase(:).cc_stimsets_all);         
   dataBase_merge.cc_stimsets_avg = cat(1,dataBase(:).cc_stimsets_avg);         
   dataBase_merge.cc_stimchans_all = cat(1,dataBase(:).cc_stimchans_all);         
   dataBase_merge.cc_stimchans_avg = cat(1,dataBase(:).cc_stimchans_avg);         
   dataBase_merge.stimpnames_all = cat(2,dataBase(:).stimpnames_all); 
   dataBase_merge.stimpnames_avg = cat(2,dataBase(:).stimpnames_avg);         
   dataBase_merge.cc_epoch_sorted = cat(3,dataBase(:).cc_epoch_sorted);     
   dataBase_merge.tt_epoch_sorted = cat(2,dataBase(:).tt_epoch_sorted);         
   dataBase_merge.cc_epoch_sorted_avg = cat(2,dataBase(:).cc_epoch_sorted_avg);         
   dataBase_merge.cc_epoch_sorted_select = cat(2,dataBase(:).cc_epoch_sorted_select);         
   dataBase_merge.stimpnames = cat(2,dataBase(:).stimpnames);  
end