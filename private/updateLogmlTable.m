function updateLogmlTable(pops)
% Updates global variables LOGML_TABLE, npops*1 array, logml values for
% each population given in "pops"
% After the updates, the values are based on the current values of the
% global variables COUNTS, SUM_COUNTS
% Lu Cheng, 11.01.2012

% transform COUNTS and SUM_COUNTS to double
% Lu Cheng
% 18.04.2012

global COUNTS;
global SUM_COUNTS; 
global LOGML_TABLE;
global PRIOR_PAR;

term1 = 0-gammaln(1+double(SUM_COUNTS(:,pops)));
term2 = gammaln(double(COUNTS(:,:,pops))+PRIOR_PAR(:,:,pops))-gammaln(PRIOR_PAR(:,:,pops));
term2(PRIOR_PAR(:,:,pops)==0) = 0;
term2 = reshape(sum(term2,1),size(term1));

LOGML_TABLE(pops) = sum(term1 + term2,1);

%----------------------------------------------------------------------