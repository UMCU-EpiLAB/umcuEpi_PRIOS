function localDataPath = personalDataPath_pros(varargin)

% function that contains local data path, is ignored in .gitignore
%localDataPath.elec_input = '/home/sifra/db/2_Chronic_ECoG/Metadata/Electrodes/';
localDataPath.elec_input = '/home/sifra/Desktop/db/respect-leijten/2_Chronic_ECoG/Metadata/Electrodes/';
localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/PRIOSstudy/'; 
localDataPath.dataPath = '/Fridge/KNF/PRIOSstudy/';
localDataPath.CCEP_allpat = '/Fridge/users/sifra/derivatives/CCEP/CCEP_files_allPat/' ;

% set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/'))     
addpath('/home/sifra/git_repositories/fieldtrip')
ft_defaults

% Remove paths (remove retrospective folder because filenames are matching)
warning('off','all');
rmpath(genpath('/home/sifra/git_repositories/CCEP_NMM_SB/CCEP/Retrospective_analysis/'))
warning('on','all');

end
