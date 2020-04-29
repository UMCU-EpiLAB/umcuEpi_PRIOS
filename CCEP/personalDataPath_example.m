
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

localDataPath.input = 'bla/CCEP/';
localDataPath.output = '/bla/users/dorien/';

% % set paths
fieldtrip_folder  = '/bla/fieldtrip/';
% copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/bla/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults

end
