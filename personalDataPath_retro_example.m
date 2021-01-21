function localDataPath = personalDataPath_retro_example(varargin)

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'sub_labels'))
                % RESPect database
                localDataPath.CCEPpath = '/folder/to/derivatives/CCEP/';
                localDataPath.dataPath = '/folder/to/bidsfiles/chronic_ECoG/';
                localDataPath.CCEP_allpat = '/folder/to/write/bidsderivatives/derivatives/CCEP/CCEP_files_allPat/' ;
        end
    end
end

localDataPath.elec_input = '/folder/to/electrodes/excel/Metadata/Electrodes';

% % set paths
addpath(genpath('/folder/to/gitrepository/eeglab/'))

fieldtrip_folder  = '/folder/to/gitrepository/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/folder/to/gitrepository/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults

% Remove paths (remove prospective folder because filenames are matching)
warning('off','all');
rmpath(genpath('/folder/to/gitrepository/CCEP_NMM/CCEP/Prospective_analysis/'))
rmpath(genpath('/folder/to/gitrepository/CCEP_NMM/NMM/'))
warning('on','all');


end
