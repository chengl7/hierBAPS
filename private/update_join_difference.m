function update_join_difference(data,logml)
% update JOIN_DIFFERENCE
% Lu Cheng, 11.01.2012

global COUNTS;      global SUM_COUNTS;
global PARTITION;
global JOIN_DIFFERENCE;

npops = size(COUNTS,3);

for i1 = 1:npops-1
    indsToBeMoved = find(PARTITION==i1);
    if isempty(indsToBeMoved)
        % Cluster i1 is empty
        JOIN_DIFFERENCE(i1,(i1+1):npops) = 0;
        JOIN_DIFFERENCE((i1+1):npops,i1) = 0;
    else
        [diffCounts diffSumCounts] = computeDiffInCounts(indsToBeMoved, data);
        
        unknown_pops = find(isnan(JOIN_DIFFERENCE(i1,(i1+1):end)));
        unknown_pops = unknown_pops+i1;
        
        COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffCounts;
        SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) - diffSumCounts;
        
        PARTITION(indsToBeMoved) = -1;
        updateLogmlTable(i1);
        
        for i2 = unknown_pops
            COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffCounts;
            SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) + diffSumCounts;        
            PARTITION(indsToBeMoved) = i2;
            
            updateLogmlTable(i2);
            logml_new = computeTotalLogml();
                
            JOIN_DIFFERENCE(i1,i2) = logml_new-logml;
            JOIN_DIFFERENCE(i2,i1) = logml_new-logml;
     
            COUNTS(:,:,i2) = COUNTS(:,:,i2) - diffCounts;
            SUM_COUNTS(:,i2) = SUM_COUNTS(:,i2) - diffSumCounts;
            PARTITION(indsToBeMoved) = i1;
            
            updateLogmlTable(i2);
        end
        
        COUNTS(:,:,i1) = COUNTS(:,:,i1) + diffCounts;
        SUM_COUNTS(:,i1) = SUM_COUNTS(:,i1) + diffSumCounts;
        
        PARTITION(indsToBeMoved) = i1;
        updateLogmlTable(i1);
    end
end