function localDataPath = personalDataPath_pros_example(varargin)

% function that contains local data path, is ignored in .gitignore
localDataPath.elec_input = '/folder/to/electrodes/excels/Metadata/Electrodes/';
localDataPath.CCEPpath = '/folder/to/location/to/save/derivatives/PRIOSstudy/'; 
localDataPath.dataPath = '/folder/with/bidsdata/PRIOSstudy/';
localDataPath.CCEP_allpat = '/folder/to/location/to/save/derivatives/CCEP_files_allPat/' ;
localDataPath.CCEP_interObVar = '/folder/to/location/with/interobserver/data/';

% set paths
addpath(genpath('/folder/to/gitrepository/eeglab/'))     
addpath('/folder/to/gitrepository/fieldtrip')
ft_defaults

% Remove paths (remove retrospective folder because filenames are matching)
warning('off','all');
rmpath(genpath('/folder/to/gitrepository/CCEP_NMM/CCEP/Retrospective_analysis/'))
rmpath(genpath('/folder/to/gitrepository/fieldtrip/external/signal')) % to avoid usage of fieldtrip butter and filtfilt instead matlab functions
warning('on','all');

end
