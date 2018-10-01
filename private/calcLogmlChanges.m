function changes = calcLogmlChanges(inds, data, logml)
% compute the logml change if the given inds are moved to another cluster
% the input inds are supposed to come from the same cluster
% changes is a npops*1 vector
% Lu Cheng, 11.01.2012

global COUNTS;
global SUM_COUNTS; 
global PARTITION;

npops = size(COUNTS,3);
changes = zeros(npops,1);
indsToBeMoved = inds;

if isempty(indsToBeMoved), return, end

i1 = PARTITION(indsToBeMoved(1));
[diffCounts diffSumCounts]= computeDiffInCounts(inds, data);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffCounts;
SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) - diffSumCounts;

% PARTITION(inds) = -1;
updateLogmlTable(i1);

for i2 = 1:npops
    if i2 ~= i1
        COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffCounts;
        SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) + diffSumCounts;
        
%         PARTITION(inds) = i2;
        updateLogmlTable(i2);
        logml_new = computeTotalLogml();
        changes(i2) = logml_new - logml;
        
        COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffCounts;
        SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) - diffSumCounts;
%         PARTITION(inds) = -1;
        updateLogmlTable(i2);
    end
end

COUNTS(:,:,i1) = COUNTS(:,:,i1) + diffCounts;
SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) + diffSumCounts;

% PARTITION(inds) = i1;
updateLogmlTable(i1);


%---------------------------------------------------------------------