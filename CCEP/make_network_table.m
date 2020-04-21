%% MAKE SURE ONLY TWO RUNS ARE COMPARED!!
% Otherwise the agreement is not determined correctly 

% This script is preceded by the pipeline_CCEP.m in which the ERs are 
% detected and saved for the 2 stims or 10 stims networks. 
%% set paths
clear 

%myDataPath = setLocalDataPath(cfg);
myDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/' ;
myDataPath.dataPath = '/Fridge/chronic_ECoG/';

% set paths
addpath(genpath('/home/sifra/git_repositories/CCEP_NMM_SB/CCEP/'));
addpath(genpath('/home/sifra/git_repositories/eeglab/')) ;    
addpath('/home/sifra/git_repositories/fieldtrip');
ft_defaults

localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; % /Fridge/users/sifra/derivatives/CCEP
localDataPath.dataPath = '/Fridge/chronic_ECoG/';

%% pre-allocation
config_CCEP

%% Load Multiple runs
files = dir(fullfile(localDataPath.CCEPpath,cfg.sub_labels{1}, 'ses-*' ,'run-*',...
    [cfg.sub_labels{1} '_ses-*_' cfg.task_label '_*'  '_CCEP_*.mat']));

if length(files) == 0
    fprintf('WARNING: No runs are found');
else
    for i = 1:size(files,1)
    names = fullfile(files(i).folder, files(i).name);
    runs(i) = load(names);
    end
end

%% Compare the electrodes which show a N1 in different networks.

if length(runs) > 1
   for k = 1:length(runs)
    N1_peak(k,:) = {runs(k).ccep.n1_peak_sample};
    stimnames(k,:) = {runs(k).ccep.stimpnames(1,:)};
    
    for n = 1:(length(stimnames)-1)   
        if length(stimnames{n,1}) ~= length(stimnames{(n+1),1})
            
            % Find the extra stimulation pairs
            C = setdiff(stimnames{n+1,1}, stimnames{n,1});
            same_stimpair1 = setdiff(stimnames{n+1,1}, C);
            
            D = setdiff(stimnames{n,1}, stimnames{n+1,1});
            same_stimpair2 = setdiff(stimnames{n,1}, D);
            
            if length(same_stimpair1) ~= length(same_stimpair2)
                fprintf('WARNING: Extra stimpairs of %s are NOT removed.', cfg.sub_labels{:});
            end
                        
            % Find the extra electrodes          
            E = setdiff(runs(n+1).ccep.ch', runs(n).ccep.ch');
            same_chan1 = setdiff(runs(n+1).ccep.ch', E);
            
            F = setdiff(runs(n).ccep.ch', runs(n+1).ccep.ch');
            same_chan2 = setdiff(runs(n).ccep.ch', F);
            
            if length(same_chan1) ~= length(same_chan2)
                fprintf('WARNING: Extra channels of %s are NOT removed.', cfg.sub_labels{:});
            end
            
            for Q = 1:length(runs)
                correct_stimpairs{Q,1} = find(ismember(stimnames{Q,1}, same_stimpair2));  
                correct_chans{Q,1} = find(ismember(runs(Q).ccep.ch', same_chan2));        
                % does not matter whether same_chan1 or 2, should be the same
                
                % Make matrix with the channels and stimpairs which are
                % present in all runs
                matrix{Q,1} = runs(Q).ccep.n1_peak_sample(correct_chans{Q}, correct_stimpairs{Q});
                adjecency_matrix{Q,1} = isnan(matrix{Q,1});  % NAN = 1, value = 0           
            end
            % For two runs, they can be added and the values will be 2, 1
            % or 0. When multiple runs are compared, new values have to be
            % filled in.
                if length(runs) == 2
                    compare_mat = adjecency_matrix{1} + adjecency_matrix{2};
                    Z = sum(compare_mat(:)==2);         % Both no-ER
                    XandY = sum(compare_mat(:)==1);     % One NAN one ER
                    W = sum(compare_mat(:)==0);         % Two ERs
                    total = W + XandY + Z ;
                    check = size(compare_mat,1) .* size(compare_mat,2);
                else
                    fprintf('WARNING: The matrices of >2 runs cannot be added up');
                end   
        end
     end  
   end
else
   N1_peak = runs.ccep.n1_peak_sample;
   stimnames = runs{:}.ccep.stimpnames;    
end

% Overall, positive and negative agreement between the matrices. 
OA = (W + Z) / (W + XandY + Z)
PA = (2 * W) / (2 * W + XandY)
NA = (2 * Z) / (2 * Z + XandY)

