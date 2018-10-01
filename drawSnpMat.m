function h = drawSnpMat(inFile,varargin)
% draw a snp matrix, each row is a DNA sequence
% partition specifies labels for each DNA sequence
%      if it is a matrix, then each row corresponds to each DNA sequence
% snpMat could be either numeric or char, but '-ACGT' should be mapped to 0:4
% additional "shuffle" parameter could added to shuffle the columns for better visualization
%     columns sorted according to [#majorityBase entropyBaseDistribution]

% Lu Cheng
% 22.03.2013

load(inFile);
snpMat = c.snpData;
heds = c.heds;
clear c

[tmpPart tmpIdx] = sortrows(partition);
snpMat=snpMat(tmpIdx,:);
partition = tmpPart(:,1);
dispHeds = heds(tmpIdx);

if nargin>1 && strcmp(varargin{1},'shuffle')
    nRow = size(snpMat,1);
    counts = histc(snpMat,0:4,1);
    freqMat = counts/nRow+1e-6;

    [~,maxInd] = max(counts,[],1);
    muNegEntropy = -sum(freqMat.*log(freqMat),1);
    
    [~,colIdx] = sortrows([maxInd; muNegEntropy]');
    snpMat = snpMat(:,colIdx);
end

if ~issorted(partition)
    warning('partition is not sorted, now we sort the partition.');
    [newPart, xidx] = sort(partition);
    partition = newPart;
    snpMat = snpMat(xidx,:);
end   
[~, yLines] = unique(partition,'first');
yLines(1) = [];

% [~, yLines2] = unique(tmpPart(:,2),'first');
% yLines2(1) = [];

h = figure;
imagesc(snpMat,[0 4]);
map = [1 1 1; 0 0 1; 0.75 0.75 0; 1 0 0; 0 0.5 0];
colormap(map);
colorbar

if ~isempty(yLines)
    hold on
    xval = xlim;
    for i=1:length(yLines)
        line(xval,[yLines(i)-0.5 yLines(i)-0.5],'LineWidth',3,'Color','k'); %yLine(i) is the center of the color block
    end
    
%     for i=1:length(yLines2)
%         line(xval,[yLines2(i)-0.5 yLines2(i)-0.5],'LineWidth',1,'Color','k');
%     end

end

if nargin>1 && strcmp(varargin{1},'shuffle')
    title('columns shuffled. white-indel, blue-A, yellow-C, red-G, green-T');
else
    title('white-indel, blue-A, yellow-C, red-G, green-T');
end

fid = fopen('figInfo.txt','w+');
fprintf(fid,'SampleID for each row in the figure:\n');
fprintf(fid,'Row\t1-layer-clusterID\tSampleID\n');
fprintf('SampleID for each row in the figure:\n');
for i=1:length(dispHeds)
    fprintf('Row %3d, 1-layer-clusterID %3d, SampleID %s;\n',i,partition(i),dispHeds{i});
    fprintf(fid,'%d\t%d\t%s\n',i,partition(i),dispHeds{i});    
end
fclose(fid);