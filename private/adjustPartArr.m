function y=adjustPartArr(partArr)
% adjust partition array for output
% each column represents the partition for the corresponding layer
% Lu Cheng, 20.03.2012

[sortPartArr, idx] = sortrows(partArr);

[nSeq nCol] = size(partArr);

preLine = sortPartArr(1,:) + 1;
labels = zeros(1,nCol);
for i=1:nSeq
    partArr(i,:) = labels + (sortPartArr(i,:)~=preLine);
    preLine = sortPartArr(i,:);
    labels = partArr(i,:);
end

y = zeros(nSeq,nCol);
y(idx,:) = partArr;