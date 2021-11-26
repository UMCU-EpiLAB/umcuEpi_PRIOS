function [dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop)
% This script is used to check whether SPES-clin and SPES-prop have the
% same electrode-stimulation pair combinations. 
% First consider SPES-clin to have more stimpairs compared to SPES-prop,
% then consider SPES-prop to have more stimpairs compared to SPES-clin. 
% Finaly also check for equal electrode names.

% First check whether the number of the stimulation pairs is equal
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
           dataBase_prop.stimpnames(:,x_all) = [];
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
   dataBase_prop.stimpnames(:,x_all) = [];
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

    % Check for same channels
    if size(dataBase_clin.ch,1) > size(dataBase_prop.ch,1)                      % If channels in clinical SPES are diff from channels propofol SPES
        diff_elec = find(~ismember(dataBase_clin.ch, dataBase_prop.ch));
        dataBase_clin.ch(diff_elec,:) = [];     
        
    elseif size(dataBase_prop.ch,1) > size(dataBase_clin.ch,1)
        diff_elec = find(~ismember(dataBase_prop.ch, dataBase_clin.ch));
        dataBase_prop.ch(diff_elec,:) = [];
    else
        if ~isequal(dataBase_clin.ch, dataBase_prop.ch)
            warning('Although number of channels are equal in both SPESclin and SPESprop, the channel names are not similar!')
        end
    end
    
    % When the channels are still unequal, print warning
        if any(~ismember(dataBase_clin.ch, dataBase_prop.ch)) || any(~ismember(dataBase_prop.ch, dataBase_clin.ch))
            warning('%s still has unequal electrodes \n',dataBase_clin.sub_label)
        end       

end
