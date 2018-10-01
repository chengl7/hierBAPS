function logml = checkLogml(alnMat, partition)
% This function calculates the logml for the given partition
% Lu Cheng
% 09.08.2013

if ischar(alnMat)
    alnMat = nt2int(alnMat);
end
cntBases = histc(alnMat,nt2int('ACGT'));
snpPosition = sum(cntBases,1)~=max(cntBases,[],1);
snpData = double(alnMat(:,snpPosition));

prior = cntBases(:,snpPosition)>0;    % added 09.08.2013, Lu Cheng

global PARTITION;
global COUNTS;
global SUM_COUNTS;
global PRIOR_PAR;
global LOGML_TABLE;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;
global DATA_TYPE;

clearGlobalVars;

nPOPS = length(unique(partition));

% Initialize PARTITION, COUNTS, SUM_COUNTS, PRIOR_PAR, PARTITION
DATA_TYPE = 'uint16'; % supports up to 65535 input sequences
if length(partition)>intmax(DATA_TYPE)
    error('DATA_TYPE is %s, maxNum=%d is less than #sequences=%d\n!',...
        DATA_TYPE,intmax(DATA_TYPE),length(partition));
end

[sumCounts, counts] = initialCounts2(partition, snpData, nPOPS);

% tmp = sum(counts,3)>0;
% sumTmp = sum(tmp,1);
% tmpPrior = repmat(1./sumTmp,5,1).*tmp;
% priorPar = repmat(tmpPrior,[1 1 nPOPS]);

prior = prior./(repmat(sum(prior,1),4,1)); % added on 02.04.2014
priorPar = repmat(prior,[1 1 nPOPS]);
COUNTS = counts; SUM_COUNTS = sumCounts;

PRIOR_PAR = priorPar;
PARTITION = partition;

clear partition counts sumCounts tmp* sumTmp priorPar;

% Initialize LOGML_TABLE:
nINDS = length(PARTITION);
LOGML_TABLE = zeros(nPOPS,1);
updateLogmlTable(1:nPOPS);

REMOVAL_DIFFERENCE = zeros(nINDS,1);
REMOVAL_DIFFERENCE(:,:) = nan;
ADDITION_DIFFERENCE = zeros(nINDS,nPOPS);
ADDITION_DIFFERENCE(:,:) = nan;
JOIN_DIFFERENCE = zeros(nPOPS, nPOPS);
JOIN_DIFFERENCE(:,:) = nan;

logml = computeTotalLogml();  % uniform K prior

clearGlobalVars;