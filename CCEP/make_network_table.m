%% pre-allocation

clear 

config_CCEP

%% set paths

%myDataPath = setLocalDataPath(cfg);
myDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/' ;
myDataPath.dataPath = '/Fridge/chronic_ECoG/';

% set paths
addpath(genpath('/home/sifra/git_repositories/eeglab/')) ;    
addpath('/home/sifra/git_repositories/fieldtrip');
ft_defaults

localDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/'; % /Fridge/users/sifra/derivatives/CCEP
localDataPath.dataPath = '/Fridge/chronic_ECoG/';


%% choose CCEP file
files = dir(fullfile(localDataPath.CCEPpath,cfg.sub_labels{1}, cfg.ses_label,'run-*',...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_CCEP.mat']));

for i = 1:size(files,1)
    names = fullfile(files(i).folder, files(i).name);
    runs(i) = load(names);
end

%% Compare the electrodes which show a N1 in different networks.

if length(runs) > 1
   for k = 1:length(runs)
    N1_peak(k,:) = {runs(k).ccep.n1_peak_sample};
    stimnames(k,:) = {runs(k).ccep.stimpnames(1,:)};
    
    for n = 1:(length(stimnames)-1)
        if length(stimnames{n,1}) ~= length(stimnames{(n+1),1})
             % Find the names of the extra stimulation pairs
            diff_stimchan = setdiff(stimnames{n,1},stimnames{n+1,1});
            index = find(ismember(runs(n).ccep.stimpnames, diff_stimchan));     

           % index(n,:) = find(ismember(runs(n).ccep.stimpnames, diff_stimchan)); % deze bestaat niet voor de kleinste array

            % Find which columns found in diff_stimchan. 
            
           
           if size(runs(n).ccep.n1_peak_sample,2) > index(1)
            
               % if ~isempty (runs(n).ccep.n1_peak_sample(:,index))          % deze bestaat niet voor de kleinste array
                  for i = index(end) :-1: index(1)
                    runs(n).ccep.n1_peak_sample(:,i) = [];
                    runs(n).ccep.n1_peak_amplitude(:,i) = [];
                    fprintf('WARNING: Extra stimpairs of %s_%s are removed.', cfg.sub_labels{:}, cfg.ses_label);
                  end 
           end 
        end
     end  
   end
else
   N1_peak = runs.ccep.n1_peak_sample;
   stimnames = runs{:}.ccep.stimpnames;    
end

%%
run = 1:length(runs)
[row,col] = find(ismember(runs(run).ccep.n1_peak_sample  , runs(run+1).ccep.n1_peak_sample))

