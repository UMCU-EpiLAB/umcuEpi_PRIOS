function agreement_parameter = agreement_parameters(agreement, dataBase10, dataBase2, myDataPath)
    agreement_parameter = struct;    
    wantedAmat10 = agreement.Amat10'; % channels x stimulus pairs
    wantedAmat2 = agreement.Amat2'; % channels x stimulus pairs
%     Amat_10 = zeros(size(wantedAmat10,1)+size(wantedAmat10,2));     % Matrix with 4 quadrants, left under must be filled with the adjacency matrix            
%     Amat_2 = zeros(size(wantedAmat2,1)+size(wantedAmat2,2));        % Matrix with 4 quadrants, left under must be filled with the adjacency matrix            

    % find start and end of left under quadrant
%     rowStart = size(wantedAmat10,2)+1;      % Same for 2 stims
%     columnEnd = size(wantedAmat10,2);
% %     ElecName = dataBase10.ch;
%        
%     Amat_10((rowStart:end),(1:columnEnd)) = wantedAmat10; % channels x stimulus pairs in left under quadrant
%     Amat_2((rowStart:end),(1:columnEnd)) = wantedAmat2;
% 
%     for stimPair = size(wantedAmat10,2)+1:size(Amat_10,1)
%         ERs_stimp10((stimPair-(rowStart-1)),1) = sum(Amat_10(stimPair,:));       % total number of ERs evoked by a stimulation pair 
%         ERs_stimp2((stimPair-(rowStart-1)),1) = sum(Amat_2(stimPair,:));
%     end
% total number of ERs evoked by a stimulation pair

ERs_stimp10 = sum(wantedAmat10,2);
ERs_stimp2 = sum(wantedAmat2,2);    

% total number of ERs evoked in an electrode
ERs_elec10 = sum(wantedAmat10,1);
ERs_elec2 = sum(wantedAmat2,1);
% for elec = 1:size(wantedAmat10,2)
%     ERs_elec10(1,elec) = sum(Amat_10(:,elec));                            % total number of ERs evoked in an electrode (el)
%     ERs_elec2(1,elec) = sum(Amat_2(:,elec));
% end
      
%% Work with the electrodes in the columns and rows 
% All stimulations (10)    
    elec_mat10 = dataBase10.elec_Amat;
    G_10 = digraph(elec_mat10);
    wbc_10 = centrality(G_10,'betweenness','Cost',G_10.Edges.Weight);           % Weighted Betweenness centrality
    Indegree_10 = centrality(G_10,'indegree','Importance',G_10.Edges.Weight);   % indegree
    Outdegree_10 = centrality(G_10,'outdegree','Importance',G_10.Edges.Weight); % outdegree
        
% Only 2 stimulations, first pos and first neg direction
    elec_mat2 = dataBase2.elec_Amat;
    G_2 = digraph(elec_mat2);
    wbc_2 = centrality(G_2,'betweenness','Cost',G_2.Edges.Weight);        
    Indegree_2 = centrality(G_2,'indegree','Importance',G_2.Edges.Weight);
    Outdegree_2 = centrality(G_2,'outdegree','Importance',G_2.Edges.Weight);
    
    
%%% Normalisation %%%
    % total number of electrodes
    stimelektot = size(dataBase10.ch,1);
    % total number of stimpairs
    stimptot = size(dataBase10.stimpnames_avg,2);
    
    % Number of times an electrode is used in a stimulation pair  
    trialelek = zeros(1,stimelektot);
    for el=1:size(elec_mat10,1)
        trialelek(el) = size(find(dataBase10.stimsets_avg==el),1);
    end
    
    % total number of possible connections
    n_outtot = zeros(1,stimelektot);
    n_intot = zeros(1,stimelektot);
    for el=1:size(elec_mat10,1)
        n_outtot(el) = trialelek(el)*(stimelektot-2); % number of stimulated x total number of elecrodes - stimulated electrodes
        n_intot(el) = 2*(stimptot - trialelek(el)); % both stimulated electrodes x (total number of stimuli - number of stimulated)
    end
    
    % NETWERKMATEN   
    outdegreenorm10 = zeros(1,stimelektot);
    indegreenorm10 = zeros(1,stimelektot);
    BCnorm10 = zeros(1,stimelektot); 
    
    outdegreenorm2 = zeros(1,stimelektot);
    indegreenorm2 = zeros(1,stimelektot);
    BCnorm2 = zeros(1,stimelektot); 
    
    for el=1:size(elec_mat10,1)
        outdegreenorm10(el) = Outdegree_10(el)/n_outtot(el);
        indegreenorm10(el) = Indegree_10(el)/n_intot(el);
        BCnorm10(el) = wbc_10(el)/(n_intot(el)*n_outtot(el));
        
        outdegreenorm2(el) = Outdegree_2(el)/n_outtot(el);
        indegreenorm2(el) = Indegree_2(el)/n_intot(el);
        BCnorm2(el) = wbc_2(el)/(n_intot(el)*n_outtot(el));
    end
    
  agreement_parameter.indegreeN_10 = indegreenorm10;
  agreement_parameter.outdegreeN_10 = outdegreenorm10;
  agreement_parameter.BCN_10 = BCnorm10;
  agreement_parameter.ERs_stimp10 = ERs_stimp10;
  agreement_parameter.ERs_elec10 = ERs_elec10;
  
  agreement_parameter.indegreeN_2 = indegreenorm2;
  agreement_parameter.outdegreeN_2 = outdegreenorm2;
  agreement_parameter.BCN_2 = BCnorm2;
  agreement_parameter.ERs_stimp2 = ERs_stimp2;
  agreement_parameter.ERs_elec2 = ERs_elec2;
  
  % All variables are also saved to a excel variant
%   write2excelTables(dataBase10, myDataPath, agreement_parameter);       % database is only needed for the electrode names and the stimnames, therefore does not matter whether database 10 or 2 is taken.
  
end


