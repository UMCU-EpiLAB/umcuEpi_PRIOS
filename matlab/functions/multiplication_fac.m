function multiplication_fac(dataBase, ccep_allPat)

% Remove patients with too low interobserver agreement
skip_pat = zeros(size(dataBase,2),1);

% Remove patients that have an interobserver agreement lower than 0.6
for i = 1:size(dataBase,2)      
   if dataBase(i).ccep_clin.Ckappa < 0.6 && dataBase(i).ccep_prop.Ckappa <0.6
       skip_pat(i,:) = 1;
   end
end

dataBase(skip_pat==1) = [];
ccep_allPat.sub_labels(skip_pat==1) = [];  

measure = {'Indegree','Outdegree','BC'};

Mult_factor = zeros(size(dataBase,2), size(measure,2));

for n = 1:size(measure,2)        

    for subj = 1:size(dataBase,2)
        
        if strcmp(measure{n},'Indegree')
             M_Clin = prctile(dataBase(subj).agreement_parameter.indegreeN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.indegreeN_Prop,[25 50 75]);
        
        elseif strcmp(measure{n},'Outdegree')
             M_Clin = prctile(dataBase(subj).agreement_parameter.outdegreeN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.outdegreeN_Prop,[25 50 75]);
        
        elseif strcmp(measure{n},'BC')
             M_Clin = prctile(dataBase(subj).agreement_parameter.BCN_Clin,[25 50 75]);
             M_Prop = prctile(dataBase(subj).agreement_parameter.BCN_Prop,[25 50 75]);
        end
    
        
        Mult_factor(subj,n) = M_Prop(2)/M_Clin(2);       
    end
end

T = table(Mult_factor(:,1),Mult_factor(:,2),Mult_factor(:,3), 'VariableNames',measure,'RowNames',ccep_allPat.sub_labels);
               
display(T)

for n=1:size(measure,2)

    Mult = sum(Mult_factor(:,n)) / 6;
    fprintf('Multiplication factor of the %s of the SPES-clin and SPES-prop = %1.1f \n', measure{n}, Mult);

end

end
