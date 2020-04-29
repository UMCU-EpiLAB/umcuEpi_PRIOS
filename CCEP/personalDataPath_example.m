
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'sub_labels'))
            if contains(varargin{1}.sub_labels,'REC2Stim')
                % REC2Stim study
                localDataPath.CCEPpath = '/XX/';
                localDataPath.dataPath = '/XX/';
            elseif contains(varargin{1}.sub_labels,'RESP')
                % RESPect database
                localDataPath.CCEPpath = '/XX/CCEP/';
                localDataPath.dataPath = '/XX/CCEP/';
            end
        end
    end
end


% % set paths
addpath(genpath('/XX/BasicCode_ECoG_DvB'))
addpath(genpath('/XX/eeglab/'))
addpath(genpath('/XX/CCEP_NMM'))

fieldtrip_folder  = '/XX/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/XX/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults

end