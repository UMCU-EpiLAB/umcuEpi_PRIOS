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
        
        else  % When channels are not equal, take the SPESclin as original               
              diff_order = find(~cellfun(@isequal, dataBase(i+1).tb_channels{:,1}, dataBase(i).tb_channels{:,1}));
              diff_clin = dataBase(i).tb_channels{diff_order,1};
              diff_prop = dataBase(i+1).tb_channels{diff_order,1};
            
              str_diff_order = [repmat('%d, ',1,size(diff_order,1)-1),'%d'];
              str_diff_elec = [repmat('%s, ',1,size(diff_order,1)-1),'%s'];

               
              error(sprintf(['On row *' str_diff_order '* in tb_channels is a difference between the two runs. \nIn the first it is called ' str_diff_elec, ...
                  '\n' 'In the second it is called ' str_diff_elec], diff_order, diff_clin{:}, diff_prop{:}))

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