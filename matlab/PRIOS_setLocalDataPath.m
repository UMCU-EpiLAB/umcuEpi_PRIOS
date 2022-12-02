
function localDataPath = PRIOS_setLocalDataPath(varargin)

% function LocalDataPath = setLocalDataPath(varargin)
% Return the path to the root directory and add paths in this repo
%
% input:
%   PRIOS_personalDataPath: optional, set to 1 if adding PRIOS_personalDataPath
%
% when adding PRIOS_personalDataPath, the following function should be in the
% root of this repo:
%
% function localDataPath = PRIOS_personalDataPath()
%     'localDataPath = [/my/path/to/data];
%
% this function is ignored in .gitignore
%

if isempty(varargin)

    rootPath = which('setLocalDataPath');
    RepoPath = fileparts(rootPath);

    % add path to functions
    addpath(genpath(RepoPath));

    % add localDataPath default
    localDataPath = fullfile(RepoPath,'data');

elseif ~isempty(varargin)
    % add path to data
    if varargin{1}==1 && exist('PRIOS_personalDataPath','file')
        
            localDataPath = PRIOS_personalDataPath;

    elseif varargin{1}==1 && ~exist('PRIOS_personalDataPath','file')

        sprintf(['add PRIOS_personalDataPath function to add your localDataPath:\n'...
            '\n'...
            'function localDataPath = PRIOS_personalDataPath()\n'...
            'localDataPath.input = [/my/path/to/data];\n'...
            'localDataPath.output = [/my/path/to/output];\n'...
            '\n'...
            'this function is ignored in .gitignore'])
    end

end


