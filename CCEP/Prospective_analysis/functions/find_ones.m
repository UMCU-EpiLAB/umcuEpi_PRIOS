function LocOnes = find_ones(dataBase, myDataPath)

   %% Location of the ones and twos
for subj = 1:size(dataBase,2)
    
    clin = dataBase(subj).ccep_clin.n1_peak_amplitude_check;
    prop = dataBase(subj).ccep_prop.n1_peak_amplitude_check;
    i = 1;
    
    for stimp = 1:size(dataBase(subj).ccep_prop.stimpnames_avg,2)                          % For each stimpair
        for elec = 1:size(dataBase(subj).ccep_prop.ch,1)                                   % For each electrode
          
            % When both clinical SPES and propofol SPES show an ER
            
              if ~isnan(clin(elec, stimp)) &&  isnan(prop(elec, stimp)) 
                   LocERclin(i,1) = dataBase(subj).ccep_prop.stimpnames_avg(stimp);
                   LocERclin(i,2) = dataBase(subj).ccep_prop.ch(elec);
                   i = i+1;
                                      
              end
        end      
    end
    
    % Write to table
    ColN = {sprintf('%s Stimpairs',dataBase(subj).sub_label), sprintf('%s Electrode',dataBase(subj).sub_label)} ;           %{'Stimpairs','Electrodes'};
    LocaTwos = array2table(LocERclin,'VariableNames',ColN);    
    
    % Save table
    targetFolder = [myDataPath.CCEPpath, 'Location ERs/'];
    fileName = ['Location_ER_clin_',dataBase(subj).sub_label,'.xlsx'];
    
    writetable(LocaTwos  ,[targetFolder, fileName],'WriteRowNames',true)

    clear LocERprop
end

end

% LocTwos(i,1) = dataBase(subj).ccep_prop.stimpnames_avg(stimp);
% LocTwos(i,2) = dataBase(subj).ccep_prop.ch(elec);
% i = i+1;
