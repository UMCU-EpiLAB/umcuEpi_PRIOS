function dataBase_merge = merge_runs(dataBase)

   dataBase_merge.sub_label =  dataBase(1).sub_label;
   dataBase_merge.ses_label =  dataBase(1).ses_label;
   dataBase_merge.task_label =  dataBase(1).task_label;
   dataBase_merge.run_label =  {dataBase(1).run_label,dataBase(2).run_label};
   dataBase_merge.dataName =  dataBase(1).dataName;
   dataBase_merge.ccep_header = dataBase(1).ccep_header;
   dataBase_merge.ch =  dataBase(1).ch;
   dataBase_merge.max_stim =  dataBase(1).max_stim;
   dataBase_merge.tt =  dataBase(1).tt;

   
   % Concatenate the stimpairs in the two runs
   dataBase_merge.cc_stimsets_all = cat(1,dataBase(1).cc_stimsets_all, dataBase(2).cc_stimsets_all);         
   dataBase_merge.cc_stimsets_avg = cat(1,dataBase(1).cc_stimsets_avg, dataBase(2).cc_stimsets_avg);         
   dataBase_merge.cc_stimchans_all = cat(1,dataBase(1).cc_stimchans_all, dataBase(2).cc_stimchans_all);         
   dataBase_merge.cc_stimchans_avg = cat(1,dataBase(1).cc_stimchans_avg, dataBase(2).cc_stimchans_avg);         
   dataBase_merge.stimpnames_all = cat(2,dataBase(1).stimpnames_all, dataBase(2).stimpnames_all); 
   dataBase_merge.stimpnames_avg = cat(2,dataBase(1).stimpnames_avg, dataBase(2).stimpnames_avg);         
   dataBase_merge.cc_epoch_sorted = cat(3,dataBase(1).cc_epoch_sorted, dataBase(2).cc_epoch_sorted);     
   dataBase_merge.tt_epoch_sorted = cat(2,dataBase(1).tt_epoch_sorted, dataBase(2).tt_epoch_sorted);         
   dataBase_merge.cc_epoch_sorted_avg = cat(2,dataBase(1).cc_epoch_sorted_avg, dataBase(2).cc_epoch_sorted_avg);         
   dataBase_merge.cc_epoch_sorted_select_avg = cat(2,dataBase(1).cc_epoch_sorted_select_avg, dataBase(2).cc_epoch_sorted_select_avg);         
   dataBase_merge.stimpnames = cat(2,dataBase(1).stimpnames, dataBase(2).stimpnames);  
   dataBase_merge.tb_channels = cat(2,dataBase(1).tb_channels, dataBase(2).tb_channels);
end