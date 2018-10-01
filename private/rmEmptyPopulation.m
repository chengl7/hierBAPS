function npops = rmEmptyPopulation()
% remove empty populations from COUNTS and SUM_COUNTS
% update PARTITION
% Lu Cheng, 11.01.2012

global COUNTS;
global SUM_COUNTS;
global PARTITION;

notEmpty = find(any(SUM_COUNTS,1));

COUNTS = COUNTS(:,:,notEmpty);
SUM_COUNTS = SUM_COUNTS(:,notEmpty);

for i=1:length(notEmpty)
    apu = (PARTITION==notEmpty(i));
    PARTITION(apu)=i;
end

npops = length(notEmpty);

