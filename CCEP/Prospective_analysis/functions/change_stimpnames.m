
function ccep_prop = change_stimpnames(ccep_prop, ccep_clin)
   
Ncount = find(ismember(ccep_clin.stimpnames_all' , ccep_prop.stimpnames_all' )==0);     % if SPESprop contains more stimpairs
names = ccep_prop.stimpnames_all(Ncount);
ElecNames = ccep_prop.ch;
StimpNames = ccep_prop.stimpnames_avg;

% Test of de stimulatiepaarnaam van SPES porp voorkomt in spes CLIN
% stop een extra 0 erbij, bij channels en bij stimulatiepaar
% Test of dit het probleem oplost
% Zo niet, dan warning geven dat stimlatieparen anders zijn.

    for i=1:size(ccep_prop.stimpnames_avg,2)   
        if  sum(~ismember(ccep_clin.stimpnames_avg{i}, ccep_prop.stimpnames_avg{i})) == 0     % Test 1 is of het orrigineel hetzelfde is  
            % Do nothing because are the same
            
        else
            letter = regexp(StimpNames{i},'[a-z,A-Z]');
            n =  length(letter);                                             % Find number of letters because only the first are necessary
            number = regexp(StimpNames{i},'[1-9]');
            test2 = {[StimpNames{i}(1:(n/2)),'0',StimpNames{i}(number(1):number(2)-1),'0',StimpNames{i}(number(2))]};  % Test 2 is of een extra 0 ervoor zorgt dat het hetzelfde is
            
            
            if sum(strcmp(ccep_clin.stimpnames_avg,test2))==1                % Test dat stimpaar voorkomt in SPESclin
                ccep_prop.stimpnames_avg(i) = test2;      
            else
                %warning('%s is not found in SPES prop',ccep_clin.stimpnames_avg{i} )
            end
        end
    end
end