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
files = dir(fullfile(localDataPath.CCEPpath,cfg.sub_labels{1}, cfg.ses_label,cfg.run_label{1},...
    [cfg.sub_labels{1} '_' cfg.ses_label '_' cfg.task_label '_*'  '_CCEP.mat']));
Name = fullfile(files(1).folder, files(1).name);
load(Name);         % CCEP file is loaded

%%
stimpairs = ccep.stimpnames';
