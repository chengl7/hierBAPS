function updateGlobalVariables(inds, i2, data)
% this function moves the samples specified by "inds" to cluser i2
% then update all the global variables, "inds" are supposed to come from the
% same cluster
% Lu Cheng, 11.01.2012

global PARTITION; 
global COUNTS; 
global SUM_COUNTS;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;
global JOIN_DIFFERENCE;

i1 = PARTITION(inds(1));
PARTITION(inds)=i2;

[diffCounts diffSumCounts]= computeDiffInCounts(inds, data);

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffCounts;
SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) - diffSumCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffCounts;
SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) + diffSumCounts;

updateLogmlTable([i1 i2]);

REMOVAL_DIFFERENCE(PARTITION==i1) = nan;
REMOVAL_DIFFERENCE(PARTITION==i2) = nan;
ADDITION_DIFFERENCE(:,[i1 i2]) = nan;

JOIN_DIFFERENCE(:,i2) = nan;
JOIN_DIFFERENCE(i2,:) = nan;

if ~any(PARTITION==i1)
    % i1 became empty
    JOIN_DIFFERENCE(:,i1) = 0;
    JOIN_DIFFERENCE(i1,:) = 0;
    JOIN_DIFFERENCE(i1,i1) = nan;
else
    JOIN_DIFFERENCE(:,i1) = nan;
    JOIN_DIFFERENCE(i1,:) = nan;
end