function update_difference_tables(ind, data, logml)
% update ADDITION_DIFFERENCE and REMOVAL_DIFFERENCE
% Lu Cheng, 11.01.2012

global COUNTS;      global SUM_COUNTS;
global PARTITION;
global ADDITION_DIFFERENCE;
global REMOVAL_DIFFERENCE;

rem_old = REMOVAL_DIFFERENCE;
add_old = ADDITION_DIFFERENCE;

[diffCounts diffSumCounts] = computeDiffInCounts(ind, data); 

i1 = PARTITION(ind);

if isnan(rem_old(ind))
    % Update removal difference for the individual:
    % note that we did NOT add the removed item to other clusters
    COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffCounts;
    SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) - diffSumCounts;
    PARTITION(ind) = -1;
    
    updateLogmlTable(i1);
    logml_new = computeTotalLogml();
    rem_old(ind) = logml_new-logml;
    
    COUNTS(:,:,i1) = COUNTS(:,:,i1) + diffCounts;
    SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) + diffSumCounts;
    PARTITION(ind) = i1;
    
    updateLogmlTable(i1);
end

new_pops = isnan(add_old(ind,:));
new_pops(i1) = 0;   % Own cluster needs never be calculated.
new_pops = find(new_pops);

for i2 = new_pops(:)'
    % Update addition differences for the individual:
    % note that we did NOT remove the item
    COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffCounts;
    SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) + diffSumCounts;    
    PARTITION(ind) = i2;
    
    updateLogmlTable(i2);
    logml_new = computeTotalLogml();
    add_old(ind,i2) = logml_new - logml;
    
    COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffCounts;
    SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) - diffSumCounts;
    PARTITION(ind) = i1;
    
    updateLogmlTable(i2);
end

REMOVAL_DIFFERENCE = rem_old;
ADDITION_DIFFERENCE = add_old;

%---------------------------------------------------------------------