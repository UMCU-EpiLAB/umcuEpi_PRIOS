% ccep_categorizeAnatomicalRegions

% temporal areas:
% G_temporal_inf, G_temporal_middle, G_temp_sup-Lateral,
% G_oc-temp_med-Parahip, G_oc-temp_lat-fusifor
roi{1} = {'37','38','34','23','21','43'};
roi_name{1} = 'temporal';
roi_temporal = [37 38 34 23 21 43];

% central areas
% G_postcentral G_precentral S_central
roi{2} = {'28','29','45','4', '3'};
roi_name{2} = 'central';
roi_central = [28 29 46 4 3];

% parietal areas:
% G_pariet_inf-Angular, G_pariet_inf-Supramar, G_parietal_sup
roi{3} = {'25','26','27', '30'};
roi_name{3} = 'parietal';
roi_parietal = [25 26 27 30];

% frontal areas:
% G_front_inf-Triangul, G_front_middle, G_front_inf-Opercular
roi{4} = {'14','15','12','16'}; % maybe add 16: G_front_sup
roi_name{4} = 'frontal';
roi_frontal = [14 15 12 16]; 


% occipital areas:
roi{5} = {'19', '20', '42 '}; 
roi_name{5} = 'occipital';
roi_occipital = [19 20 42 ]; 



