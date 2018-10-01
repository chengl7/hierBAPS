function logml = computeTotalLogml
% compute the log marginal likelihood of the data
% Here we use uniform K prior

% Lu Cheng, 23.01.2012

global LOGML_TABLE;
global PARTITION;

% logml = sum(LOGML_TABLE);

nCluster = length(unique(PARTITION(PARTITION>0)));
nInds = sum(PARTITION>0);
logml = sum(LOGML_TABLE)-lnStirlingS2(nInds,nCluster);