function [dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop)
% This script is used to check whether SPES-clin and SPES-prop have the
% same electrode-stimulation pair combinations. 
% First check for equal electrode names 
% Then consider SPES-clin to have more stimpairs compared to SPES-prop,
% then consider SPES-prop to have more stimpairs compared to SPES-clin. 

%% Check for same channels
% If arrays are still not equal, check the order of the electrodes
if ~isequal(dataBase_clin.ch, dataBase_prop.ch) 
    diff_order = find(~cellfun(@isequal, dataBase_prop.ch, dataBase_clin.ch));
    diff_clin = dataBase_clin.ch(diff_order);
    diff_prop = dataBase_prop.ch(diff_order);

    str_diff = [repmat('%s, ',1,size(diff_order,1)-1),'%s'];
    str_loc = [repmat('%d, ',1,size(diff_order,1)-1),'%d'];
    
    error(sprintf(['On row *' str_loc '* in the channels array is a difference between the two runs. \nIn SPESclin it is called ' str_diff, ...
        '\n' 'In SPESprop it is called ' str_diff], diff_order, diff_clin{:}, diff_prop{:}))
end


%% Check whether the number of the stimulation pairs is equal
if length(dataBase_clin.stimpnames_all) > length(dataBase_prop.stimpnames_all)                    
    [x_all,~ ] = find(ismember(dataBase_clin.stimpnames_all' , dataBase_prop.stimpnames_all' )==0); % Find the stimulation pairs extra in SPES-clin 
    [x_avg,~] = find(ismember(dataBase_clin.stimpnames_avg' , dataBase_prop.stimpnames_avg' )==0) ;
   
    names = dataBase_clin.stimpnames_all(x_all);

    stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
    fprintf(['Stimpairs only stimulated in SPESclin and not in SPESprop: \n' stringsz '\n'],names{:})

    % Remove the stimulation pairs that are only present in SPES-clin
    dataBase_clin.cc_stimsets_all(x_all,:) = [];
    dataBase_clin.cc_stimchans_all(x_all,:) = [];
    dataBase_clin.stimpnames_all(:,x_all) = [];
    dataBase_clin.cc_epoch_sorted(:,:,x_all,:) = [];
    dataBase_clin.tt_epoch_sorted(:,x_all,:) = [];
    
    
    dataBase_clin.cc_stimsets_avg(x_avg,:) = [];
    dataBase_clin.cc_stimchans_avg(x_avg,:) = [];
    dataBase_clin.stimpnames_avg(x_avg) = [];
    dataBase_clin.cc_epoch_sorted_avg(:,x_avg,:) = [];
    dataBase_clin.cc_epoch_sorted_select_reref(:,x_avg,:,:) = [];
    dataBase_clin.cc_epoch_sorted_select_reref_avg(:,x_avg,:) = [];
    dataBase_clin.cc_epoch_sorted_select(:,x_avg,:,:) = [];
    
    % Check again if previous step was enough.
    if length(dataBase_clin.stimpnames_all) ~=  length(dataBase_prop.stimpnames_all) 
        if length(dataBase_prop.stimpnames_all) > length(dataBase_clin.stimpnames_all)
           [x_all,~] = find(ismember(dataBase_prop.stimpnames_all' , dataBase_clin.stimpnames_all' )==0);     % if SPESprop contains more stimpairs
           [x_avg,~] = find(ismember(dataBase_prop.stimpnames_avg' , dataBase_clin.stimpnames_avg' )==0);     
           
           names = dataBase_prop.stimpnames_all(x_all);
           stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
           fprintf(['Stimpairs only stimulated in SPESprop and not in SPESclin: \n' stringsz '\n'],names{:})
            
           % Remove the stimulation pairs that are only present in SPES-prop
           dataBase_prop.cc_stimsets_all(x_all,:) = [];
           dataBase_prop.cc_stimchans_all(x_all,:) = [];
           dataBase_prop.stimpnames_all(:,x_all) = [];
           dataBase_prop.cc_epoch_sorted(:,:,x_all,:) = [];
           dataBase_prop.tt_epoch_sorted(:,x_all,:) = [];

           dataBase_prop.cc_stimsets_avg(x_avg,:) = [];
           dataBase_prop.cc_stimchans_avg(x_avg,:) = [];
           dataBase_prop.stimpnames_avg(x_avg) = [];
           dataBase_prop.cc_epoch_sorted_avg(:,x_avg,:) = [];
           dataBase_prop.cc_epoch_sorted_select_reref(:,x_avg,:,:) = [];
           dataBase_prop.cc_epoch_sorted_select_reref_avg(:,x_avg,:) = [];
           dataBase_prop.cc_epoch_sorted_select(:,x_avg,:,:) = [];
            
        end
    end
        % When the stimulation pairs are still unequal, print warning
        if sum(ismember(dataBase_clin.stimpnames_all, dataBase_prop.stimpnames_all)) ~= size(dataBase_clin.stimpnames_all,2)  || sum(ismember(dataBase_clin.stimpnames_all, dataBase_prop.stimpnames_all)) ~= size(dataBase_prop.stimpnames_all,2)
            warning('%s still has unequal stimulation pairs.. \n',dataBase_clin.sub_label)
        end
        
elseif length(dataBase_prop.stimpnames_all) > length(dataBase_clin.stimpnames_all)
   [x_all,~] = find(ismember(dataBase_prop.stimpnames_all' , dataBase_clin.stimpnames_all' )==0);     % if SPESprop contains more stimpairs
   [x_avg,~] = find(ismember(dataBase_prop.stimpnames_avg' , dataBase_clin.stimpnames_avg' )==0);     % if SPESprop contains more stimpairs
 
   names = dataBase_prop.stimpnames_all(x_all);
   stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
   sprintf(['Stimpairs only stimulated in SPESprop and not in SPESclin: \n' stringsz '\n'],names{:})
   
   % Remove the stimulation pairs that are only present in SPES-prop    
   dataBase_prop.cc_stimsets_all(x_all,:) = [];
   dataBase_prop.cc_stimchans_all(x_all,:) = [];
   dataBase_prop.stimpnames_all(:,x_all) = [];
   dataBase_prop.cc_epoch_sorted(:,:,x_all,:) = [];
   dataBase_prop.tt_epoch_sorted(:,x_all,:) = [];
   
   dataBase_prop.cc_stimsets_avg(x_avg,:) = [];
   dataBase_prop.cc_stimchans_avg(x_avg,:) = [];
   dataBase_prop.stimpnames_avg(x_avg) = [];
   dataBase_prop.cc_epoch_sorted_avg(:,x_avg,:) = [];
   dataBase_prop.cc_epoch_sorted_select_reref(:,x_avg,:,:) = [];
   dataBase_prop.cc_epoch_sorted_select_reref_avg(:,x_avg,:) = [];
   dataBase_prop.cc_epoch_sorted_select(:,x_avg,:,:) = [];
    
    % Check if the same number of stimpairs are used now.
    if length(dataBase_clin.stimpnames_all) ==  length(dataBase_prop.stimpnames_all) 
        if length(dataBase_clin.stimpnames_all) > length(dataBase_prop.stimpnames_all)
             [x_all,~ ] = find(ismember(dataBase_clin.stimpnames_all' , dataBase_prop.stimpnames_all' )==0); 
            [x_avg,~] = find(ismember(dataBase_clin.stimpnames_avg' , dataBase_prop.stimpnames_avg' )==0) ;

            names = dataBase_clin.stimpnames_all(x_all);

            stringsz = [repmat('%s, ',1,size(names,2)-1),'%s'];
            sprintf(['Stimpairs only stimulated in SPESclin and not in SPESprop: \n' stringsz '\n'],names{:})
            
            % Remove the stimulation pairs that are only present in SPES-clin
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
            dataBase_clin.cc_epoch_sorted_select_reref(:,x_avg,:,:) = [];
            dataBase_clin.cc_epoch_sorted_select_reref_avg(:,x_avg,:) = [];
            dataBase_clin.cc_epoch_sorted_select(:,x_avg,:,:) = [];
    
        end
    end

        % When the stimulation pairs are still unequal, print warning
        if sum(ismember(dataBase_clin.stimpnames_all, dataBase_prop.stimpnames_all)) ~= size(dataBase_clin.stimpnames_all,2)  || sum(ismember(dataBase_clin.stimpnames_all, dataBase_prop.stimpnames_all)) ~= size(dataBase_prop.stimpnames_all,2)
            warning('%s still has unequal stimulation pairs.. \n',dataBase_clin.sub_label)
        end
        
end

end
