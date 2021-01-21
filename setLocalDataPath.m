
function localDataPath = setLocalDataPath(varargin)

% function LocalDataPath = setLocalDataPath(varargin)
% Return the path to the root CCEP  directory and add paths in this repo
%
% input:
%   personalDataPath: optional, set to 1 if adding personalDataPath
%
% when adding personalDataPath, the following function should be in the
% root of this repo:
%
% function localDataPath = personalDataPath()
%     'localDataPath = [/my/path/to/data];
%
% this function is ignored in .gitignore
%
% dhermes, 2020, Multimodal Neuroimaging Lab

if isempty(varargin)
    
    rootPath = which('setLocalDataPath');
    ccepRepoPath = fileparts(rootPath);
    
    % add path to functions
    addpath(genpath(ccepRepoPath));
    
    % add localDataPath default
    localDataPath = fullfile(ccepRepoPath,'data');
    
elseif ~isempty(varargin)
    % add path to data
    if isstruct(varargin{1})
        
        if strcmp(varargin{1}.mode,'retro')
                       
            rootPath = which('setLocalDataPath');
            RepoPath = fileparts(rootPath);
                        
            % add path to functions
            addpath(genpath([RepoPath,'/CCEP/Retrospective_analysis']));
            
            localDataPath = personalDataPath_retro(varargin{1});
            
        elseif strcmp(varargin{1}.mode,'pros')
           
            rootPath = which('setLocalDataPath');
            RepoPath = fileparts(rootPath);
            
            % add path to functions
            addpath(genpath(RepoPath));
            
            % add path to functions
            addpath(genpath([RepoPath,'/CCEP/Prospective_analysis']));

            localDataPath = personalDataPath_pros(varargin{1});
            
        elseif  strcmp(varargin{1}.mode,'NMM')
            
            rootPath = which('setLocalDataPath');
            RepoPath = fileparts(rootPath);
            
            % add path to functions
            addpath(genpath(RepoPath));

        end
        
    else
        if varargin{1}==1 && exist('personalDataPath','file')
            
            if strcmp(varargin{1}.mode,'retro')
                localDataPath = personalDataPath_retro(varargin{1});
                
            elseif strcmp(varargin{1}.mode,'pros')
                localDataPath = personalDataPath_pros(varargin{1});
                
            end
            
        elseif varargin{1}==1 && ~exist('personalDataPath','file')
            
            sprintf(['add personalDataPath function to add your localDataPath:\n'...
                '\n'...
                'function localDataPath = personalDataPath()\n'...
                'localDataPath.input = [/my/path/to/data];\n'...
                'localDataPath.output = [/my/path/to/output];\n'...
                '\n'...
                'this function is ignored in .gitignore'])
            return
        end
    end 
end

return

