% This script is used to check whether SPES-clin and SPES-prop have the
% same electrode-stimulation pair combinations. 
% First check for equal electrode names, then compare the stimulus pairs
% between SPES-clin and SPES-prop.

% INPUTS:
% - dataBase_clin and dataBase_prop
%   both structs with the same fields:
%   - ch
%   cell[channels x 1] containing the channel names 
%   - cc_stimchans
%   cell[stim pairs x 2] containing all the electrode names that are
%   stimulated (C1-C2 and C2-C1 are combined)
%   - cc_stimsets
%   matrix[stim pairs x 2] containing all the electrode numbers that are
%   stimulated (C1-C2 and C2-C1 are combined)
%   - cc_epoch_sorted
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined)
%   - cc_epoch_sorted_avg
%   matrix[channels x stimulus pairs x samples] containing averaged
%   responses to stimulation (averaged from cc_epoch_sorted)
%   - cc_epoch_sorted_reref
%   matrix[channels x trials x stimulus pairs x samples] containing
%   responses to stimulation of all stimulus pairs (C01-C02 and C02-C01
%   combined) --> re-referenced
%   - cc_epoch_sorted_reref_avg
%   matrix[channels x stimulus pairs x samples] containing averaged
%   responses to stimulation (averaged from cc_epoch_sorted) -->
%   re-referenced
%   - tt_epoch_sorted
%   matrix[trials x stimulus pairs x samples] containing epoched time
%   points

% OUTPUTS:
% - dataBase:
%   the fields in this struct, mentioned in INPUTS, are adapted so that
%   only the stimulus pairs that are stimulated in both SPES-clin and
%   SPES-prop remain.

function [dataBase_clin, dataBase_prop] = similar_stimpairs(dataBase_clin, dataBase_prop)

%% Check for same channels
if ~isequal(dataBase_clin.ch, dataBase_prop.ch) 
    diff_order = find(~cellfun(@isequal, dataBase_prop.ch, dataBase_clin.ch));
    diff_clin = dataBase_clin.ch(diff_order);
    diff_prop = dataBase_prop.ch(diff_order);

    error(['For example: On row %d, %d in the channels array, a ',...
        'difference is found between the two runs. \nIn SPESclin, ',...
        'it is called %s, %s, in SPESprop %s, %s. \n',...
        'CHECK FILES FOR OTHER DIFFERENCES!\n'], ...
        diff_order(1), diff_order(2), diff_clin{1}, ...
        diff_clin{2}, diff_prop{1}, diff_prop{2})

end

%% Check whether the same stimulation pairs are stimulated, and otherwise
% delete stimulation pairs that are only stimulated in either propofol or
% clinical SPES

stimPropPresent = zeros(size(dataBase_prop.cc_stimchans,1),1);
stimClinPresent = zeros(size(dataBase_clin.cc_stimchans,1),1);

for iStim = 1:size(dataBase_clin.cc_stimchans,1)
    iStimProp = find(strcmpi(dataBase_prop.cc_stimchans(:,1),dataBase_clin.cc_stimchans{iStim,1}) & ...
        strcmpi(dataBase_prop.cc_stimchans(:,2),dataBase_clin.cc_stimchans{iStim,2}));
    
    if ~isempty(iStimProp)
        stimClinPresent(iStim) = 1;
        stimPropPresent(iStimProp) = 1;
    end

end

% Remove the stimulation pairs that are only present in SPES-clin
dataBase_clin.cc_stimsets(stimClinPresent == 0,:) = [];
dataBase_clin.cc_stimchans(stimClinPresent == 0,:) = [];
dataBase_clin.cc_epoch_sorted(:,stimClinPresent == 0,:,:) = [];
dataBase_clin.cc_epoch_sorted_avg(:,stimClinPresent == 0,:) = [];
dataBase_clin.cc_epoch_sorted_reref(:,stimClinPresent == 0,:,:) = [];
dataBase_clin.cc_epoch_sorted_reref_avg(:,stimClinPresent == 0,:) = [];
dataBase_clin.tt_epoch_sorted(:,stimClinPresent == 0,:) = [];

% Remove the stimulation pairs that are only present in SPES-prop
dataBase_prop.cc_stimsets(stimPropPresent == 0,:) = [];
dataBase_prop.cc_stimchans(stimPropPresent == 0,:) = [];
dataBase_prop.cc_epoch_sorted(:,stimPropPresent == 0,:,:) = [];
dataBase_prop.cc_epoch_sorted_avg(:,stimPropPresent == 0,:) = [];
dataBase_prop.cc_epoch_sorted_reref(:,stimPropPresent == 0,:,:) = [];
dataBase_prop.cc_epoch_sorted_reref_avg(:,stimPropPresent == 0,:) = [];
dataBase_prop.tt_epoch_sorted(:,stimPropPresent == 0,:) = [];

% When the stimulation pairs are still unequal, print warning
if ~isequal(dataBase_clin.cc_stimsets,dataBase_prop.cc_stimsets)
    warning('%s still has unequal stimulation pairs.. \n',dataBase_clin.sub_label)
end

end
