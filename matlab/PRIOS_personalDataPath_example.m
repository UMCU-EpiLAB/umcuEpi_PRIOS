function localDataPath = PRIOS_personalDataPath_example(varargin)
% This function contains local data path, and is ignored in .gitignore

% \\\
% ATTENTION: THIS FUNCTION IS AN EXAMPLE!
% \\\

% You should make your own PRIOS_personalDataPath.m where you fill in the
% correct repositories. This PRIOS_personalDataPath.m is ignored in
% .gitignore and should NEVER be visible online! Only
% PRIOS_personalDataPath_example.m should be visible online. 

% load BIDS data from:
localDataPath.dataPath = '/folder/to/BIDSdata/'; 
% save processed data here:
localDataPath.CCEPpath = '/folder/to/save/derivatives/of/priosstudy/'; 
% save figures here:
localDataPath.Figures = '/folder/to/save/figures/of/results/';

% set paths
fieldtrip_folder  = '/folder/to/gitrepository/fieldtrip/';
fieldtrip_private = '/folder/to/gitrepository/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults

end
