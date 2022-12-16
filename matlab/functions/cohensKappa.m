% calculate unweighted cohen's kappa

% INPUT:
% - R1_obs
%   vector[observations x 1] with 0 (no observation) and 1 (observation) of
%   observer 1
% - R2_obs
%   vector[observations x 1] with 0 (no observation) and 1 (observation) of
%   observer 2

% OUTPUT:
% kappa: scalar with the calculated unweighted Cohen's kappa

function kappa = cohensKappa(R1_obs, R2_obs)

C = confusionmat(R1_obs, R2_obs);          %disp(C);            % Convert to confusion matrix
n = sum(C(:));                                                  % get total N
C = C./n;                                                       % Convert confusion matrix counts to proportion of n
r = sum(C,2);                                                   % row sum
s = sum(C);                                                     % column sum
expected = r*s;                                                 % expected proportion for random agree
po = sum(diag(C));                                              % Observed proportion correct
pe = sum(diag(expected));                                       % Proportion correct expected
kappa = (po-pe)/(1-pe);                                         % Cohen's kappa

end