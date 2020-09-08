function localDataPath = personalDataPath_retro(varargin)

% function that contains local data path, is ignored in .gitignore
localDataPath.elec_input = '/home/sifra/db/2_Chronic_ECoG/Metadata/Electrodes/';
localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; 
localDataPath.dataPath = '/Fridge/KNF/chronic_ECoG/';
localDataPath.CCEP_allpat = '/Fridge/users/sifra/derivatives/CCEP/CCEP_files_allPat/' ;

% set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/'))     
addpath('/home/sifra/git_repositories/fieldtrip')
ft_defaults

% Remove paths (remove prospective folder because filenames are matching)
warning('off','all');
rmpath(genpath('/home/sifra/git_repositories/CCEP_NMM_SB/CCEP/Prospective_analysis/'))
warning('on','all');

end