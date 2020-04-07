
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin(1))
    if isstruct(varargin(1))
        if sum(contains(fieldnames(varargin(1)),'sub_labels'))
            if contains(varargin(1).sub_labels,'RESP')
                % RESPect database
                %localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; 
                %localDataPath.dataPath = '/Fridge/chronic_ECoG/';
                CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; 
                dataPath = '/Fridge/chronic_ECoG/';
            end
        end
    end
end


% % set paths
addpath(genpath('/home/sifra/git_repositories/CCEP_NMM_SB/CCEP/Dorien new'))
addpath(genpath('/home/sifra/git_repositories/eeglab/'))     
addpath('/home/sifra/git_repositories/fieldtrip')
ft_defaults

end