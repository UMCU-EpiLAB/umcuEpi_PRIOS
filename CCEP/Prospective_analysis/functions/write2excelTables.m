function write2excelTables(dataBase, myDataPath, agreement_parameter)

% Sorted table with the number of ERs per stimulation pair

%%% HIER VERDER MET 10 OMZETTTEN NAAR CLIN EN 2 OMZETTEN NAAR PROP

    ColN = {'number of ERs'};
    rowNames = dataBase.stimpnames_avg' ;
    rank_stimp10 = array2table(agreement_parameter.ERs_stimp10,'RowNames',rowNames,'VariableNames',ColN);   % Number of outgoing stimulations (ERs evoked per stimpair)
    rank_stimp10 = sortrows(rank_stimp10,1,'descend');
    [~, ~, ic] = unique(rank_stimp10);
    Ranked = num2cell(max(ic)-ic+1) ;
    rank_stimp10(:,2) = Ranked;
    rank_stimp10.Properties.VariableNames{2} = 'Assigned loc';   
    
    rank_stimp2 = array2table(agreement_parameter.ERs_stimp2,'RowNames',rowNames,'VariableNames',ColN);   
    rank_stimp2 = sortrows(rank_stimp2,1,'descend');
    [~, ~, ic] = unique(rank_stimp2);
    Ranked = num2cell(max(ic)-ic+1) ;
    rank_stimp2(:,2) = Ranked;
    rank_stimp2.Properties.VariableNames{2} = 'Assigned loc'; 
    
% Sorted table with the number of ERs per electrode
    ERs_elec10 = agreement_parameter.ERs_elec10';
    ERs_elec2 = agreement_parameter.ERs_elec2';
    RowNames = dataBase.ch; 
    rank_elec10 = array2table(ERs_elec10,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    rank_elec10 = sortrows(rank_elec10,1,'descend');
    rank_elec2 = array2table(ERs_elec2,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    rank_elec2 = sortrows(rank_elec2,1,'descend');
    
    
% Sorted table with the normalised indegree
    indegreeN_10 = agreement_parameter.indegreeN_10' ;
    indegreeN_2 = agreement_parameter.indegreeN_2' ;
    
    indeg10 = array2table(indegreeN_10 ,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    indeg10 = sortrows(indeg10,1,'descend');
    indeg2 = array2table(indegreeN_2,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    indeg2 = sortrows(indeg2,1,'descend');
    
    
% Sorted table with the normalised outdegree
    outdegreeN_10 = agreement_parameter.outdegreeN_10' ;
    outdegreeN_2 = agreement_parameter.outdegreeN_2' ;
    
    outdeg10 = array2table(outdegreeN_10 ,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    outdeg10 = sortrows(outdeg10,1,'descend');
    outdeg2 = array2table(outdegreeN_2,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    outdeg2 = sortrows(outdeg2,1,'descend');
    
 % Sorted table with the normalised BC
    BCN_10 = agreement_parameter.BCN_all_10' ;
    BCN_2 = agreement_parameter.BCN_2' ;
    
    BC10 = array2table(BCN_10 ,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    BC10 = sortrows(BC10,1,'descend');
    BC2 = array2table(BCN_2,'RowNames',RowNames,'VariableNames',ColN);   % Indegree
    BC2 = sortrows(BC2,1,'descend');
   
    sub_name = [extractBetween(dataBase.dataName,'ieeg/','_ses')];
    
   
    
    % Write all to one table and save in folders
    targetFolder = [myDataPath.CCEPpath, 'Agreement_par_tables/'];
    fileName = ['Agreement_table_',sub_name{1},'.xlsx'];
    
    sheet = 'rank_stimp10';
    writetable(rank_stimp10  ,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'rank_stimp2';
    writetable(rank_stimp2 ,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet =  'rank_elec10';
    writetable(rank_elec10,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'rank_elec2';
    writetable(rank_elec2 ,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'Nindeg10';
    writetable(indeg10 ,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'Nindeg2';
    writetable(indeg2  ,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet =  'Noutdeg10';
    writetable(outdeg10,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'Noutdeg2';
    writetable(outdeg2,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet =  'BC10';
    writetable(BC10,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

    sheet = 'BC2';
    writetable(BC2,[targetFolder, fileName],'sheet',sheet,'WriteRowNames',true)

end
