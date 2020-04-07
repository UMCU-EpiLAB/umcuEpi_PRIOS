%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin)
    if isstruct(varargin)
        if sum(contains(fieldnames(varargin),'sub_labels'))
            localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; % /Fridge/users/sifra/derivatives/CCEP
            localDataPath.dataPath = '/Fridge/chronic_ECoG/';
        end
    end
end

% set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/'))     
addpath('/home/sifra/git_repositories/fieldtrip')
ft_defaults

disp('This worked')

end