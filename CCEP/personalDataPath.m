function localDataPath = personalDataPath(cfg)

% function that contains local data path, is ignored in .gitignore

if ~isempty(cfg(1))
    if isstruct(cfg(1))
        if sum(contains(fieldnames(cfg(1)),'sub_labels'))
            if contains(cfg(1).sub_labels,'RESP')
                localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; 
                localDataPath.dataPath = '/Fridge/chronic_ECoG/';
            end
        end
    end
end


% % set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/'))     
addpath('/home/sifra/git_repositories/fieldtrip')
ft_defaults
end