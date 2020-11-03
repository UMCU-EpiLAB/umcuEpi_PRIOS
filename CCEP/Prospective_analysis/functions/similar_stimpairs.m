function [dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop)

if length(dataBase_clin.stimpnames_all) > length(dataBase_prop.stimpnames_all)                    % if SPESclin contains more stimpairs
    [x_all,~ ] = find(ismember(dataBase_clin.stimpnames_all' , dataBase_prop.stimpnames_all' )==0); 
    [x_avg,~] = find(ismember(dataBase_clin.stimpnames_avg' , dataBase_prop.stimpnames_avg' )==0) ;
   
    names = dataBase_clin.stimpnames_all(x_all);

    stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
    sprintf(['Stimpairs only stimulated in SPESclin and not in SPESprop: \n' stringsz '\n'],names{:})

    dataBase_clin.cc_stimsets_all(x_all,:) = [];
    dataBase_clin.cc_stimchans_all(x_all,:) = [];
    dataBase_clin.stimpnames_all(:,x_all) = [];
    dataBase_clin.stimpnames(:,x_all) = [];
    dataBase_clin.cc_epoch_sorted(:,:,x_all,:) = [];
    dataBase_clin.tt_epoch_sorted(:,x_all,:) = [];
    
    dataBase_clin.cc_stimsets_avg(x_avg,:) = [];
    dataBase_clin.cc_stimchans_avg(x_avg,:) = [];
    dataBase_clin.stimpnames_avg(x_avg) = [];
    dataBase_clin.cc_epoch_sorted_avg(:,x_avg,:) = [];
    dataBase_clin.cc_epoch_sorted_select_avg(:,x_avg,:,:) = [];
    
    % Check if the same number of stimpairs are used now.
    if length(dataBase_clin.stimpnames_all) ~=  length(dataBase_prop.stimpnames_all) 
        if length(dataBase_prop.stimpnames_all) > length(dataBase_clin.stimpnames_all)
           [x_all,~] = find(ismember(dataBase_prop.stimpnames_all' , dataBase_clin.stimpnames_all' )==0);     % if SPESprop contains more stimpairs
           [x_avg,~] = find(ismember(dataBase_prop.stimpnames_avg' , dataBase_clin.stimpnames_avg' )==0);     % if SPESprop contains more stimpairs

           names = dataBase_prop.stimpnames_all(x_all);
           stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
           sprintf(['Stimpairs only stimulated in SPESprop and not in SPESclin: \n' stringsz '\n'],names{:})

           dataBase_prop.cc_stimsets_all(x_all,:) = [];
           dataBase_prop.cc_stimchans_all(x_all,:) = [];
           dataBase_prop.stimpnames_all(:,x_all) = [];
           dataBase_prop.stimpnames(:,x_all) = [];
           dataBase_prop.cc_epoch_sorted(:,:,x_all,:) = [];
           dataBase_prop.tt_epoch_sorted(:,x_all,:) = [];

           dataBase_prop.cc_stimsets_avg(x_avg,:) = [];
           dataBase_prop.cc_stimchans_avg(x_avg,:) = [];
           dataBase_prop.stimpnames_avg(x_avg) = [];
           dataBase_prop.cc_epoch_sorted_avg(:,x_avg,:) = [];
           dataBase_prop.cc_epoch_sorted_select_avg(:,x_avg,:,:) = [];
            
        end
    end
    
        if length(dataBase_clin.stimpnames_all) ~=  length(dataBase_prop.stimpnames_all) 
            fprintf('%s still has unequal stimulation pairs.. \n',dataBase_clin.sub_label)
        end
        
elseif length(dataBase_prop.stimpnames_all) > length(dataBase_clin.stimpnames_all)
   [x_all,~] = find(ismember(dataBase_prop.stimpnames_all' , dataBase_clin.stimpnames_all' )==0);     % if SPESprop contains more stimpairs
   [x_avg,~] = find(ismember(dataBase_prop.stimpnames_avg' , dataBase_clin.stimpnames_avg' )==0);     % if SPESprop contains more stimpairs
 
   names = dataBase_prop.stimpnames_all(x_all);
   stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
   sprintf(['Stimpairs only stimulated in SPESprop and not in SPESclin: \n' stringsz '\n'],names{:})
    
   dataBase_prop.cc_stimsets_all(x_all,:) = [];
   dataBase_prop.cc_stimchans_all(x_all,:) = [];
   dataBase_prop.stimpnames_all(:,x_all) = [];
   dataBase_prop.stimpnames(:,x_all) = [];
   dataBase_prop.cc_epoch_sorted(:,:,x_all,:) = [];
   dataBase_prop.tt_epoch_sorted(:,x_all,:) = [];
   
   dataBase_prop.cc_stimsets_avg(x_avg,:) = [];
   dataBase_prop.cc_stimchans_avg(x_avg,:) = [];
   dataBase_prop.stimpnames_avg(x_avg) = [];
   dataBase_prop.cc_epoch_sorted_avg(:,x_avg,:) = [];
   dataBase_prop.cc_epoch_sorted_select_avg(:,x_avg,:,:) = [];
    
    % Check if the same number of stimpairs are used now.
    if length(dataBase_clin.stimpnames_all) ==  length(dataBase_prop.stimpnames_all) 
        if length(dataBase_clin.stimpnames_all) > length(dataBase_prop.stimpnames_all)
             [x_all,~ ] = find(ismember(dataBase_clin.stimpnames_all' , dataBase_prop.stimpnames_all' )==0); 
            [x_avg,~] = find(ismember(dataBase_clin.stimpnames_avg' , dataBase_prop.stimpnames_avg' )==0) ;

            names = dataBase_clin.stimpnames_all(x_all);

            stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
            sprintf(['Stimpairs only stimulated in SPESclin and not in SPESprop: \n' stringsz '\n'],names{:})

            dataBase_clin.cc_stimsets_all(x_all,:) = [];
            dataBase_clin.cc_stimchans_all(x_all,:) = [];
            dataBase_clin.stimpnames_all(:,x_all) = [];
            dataBase_clin.stimpnames(:,x_all) = [];
            dataBase_clin.cc_epoch_sorted(:,:,x_all,:) = [];
            dataBase_clin.tt_epoch_sorted(:,x_all,:) = [];

            dataBase_clin.cc_stimsets_avg(x_avg,:) = [];
            dataBase_clin.cc_stimchans_avg(x_avg,:) = [];
            dataBase_clin.stimpnames_avg(x_avg) = [];
            dataBase_clin.cc_epoch_sorted_avg(:,x_avg,:) = [];
            dataBase_clin.cc_epoch_sorted_select_avg(:,x_avg,:,:) = [];
    
        end
    end

    
        if length(dataBase_clin.stimpnames_all) ~=  length(dataBase_prop.stimpnames_all) 
            fprintf('%s still has unequal stimulation pairs.. \n',dataBase_clin.sub_label)
        end
end
end