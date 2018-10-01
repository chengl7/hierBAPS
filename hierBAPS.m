function hierBAPS(inFile,maxDepth, nPops, outFile)
% This function performs hierarchical BAPS clustering to the sequence alignment
% alnMat: alignment of DNA sequences, should contain only 'A','C','G','T','-'
% maxDepth: maximum depth of hierarchical search
% nPops: maximum number of population in the data
% Lu Cheng
% 19.01.2012

% REUSE_MODE was added, we can input former result c to alnMat to continue running several layers.
% Quit the process if the first layer detect the same number of clusters as nPops
% Lu Cheng
% 18.03.2012

% adjustPartArr.m added, so the output partition are in a better format
% 20.03.2012

% updated 09.08.2013, regarding the SNP site identification and prior specification
% Lu Cheng

% updated 16.01.2015, add the '-v7.3' flag to save command
% Lu Cheng

if ischar(maxDepth)
    maxDepth = str2num(maxDepth);
end

if ischar(nPops)
    nPops = str2num(nPops);
end

varName = who('-file',inFile);

REUSE_MODE = false;

if strcmp(varName{1},'alnMat')
    load(inFile);
    if size(alnMat,1)<3
        fprintf('less than 3 sequences. Quit!\n');
        return;
    end

    if ~all(ismember(unique(alnMat(:)),'ACGT-'))
        error('Unkown character in the alignment: %s. Must be ACGT-.\n',unique(alnMat(:))');
    end
elseif strcmp(varName{1},'c')
    REUSE_MODE = true;
    load(inFile);
    if c.maxDepth>=maxDepth
        error('Preprocessed data maxDepth=%d >= current maxDepth=%d. Set a larger maxDepth!\b',c.maxDepth,maxDepth);
    end
    tmpArr = zeros(c.nOrigSeq,maxDepth);
    tmpArr(:,1:c.depth)=c.partArr;
    c.partArr = tmpArr;
    c.maxDepth = maxDepth;
    clear alnMat tmpArr
else
    error('Wrong input file %s.\n',inFile);
end

if ~isnumeric(maxDepth) || maxDepth<=0
    error('maxDepth must be a number higher than 0. maxDepth=%d\n',maxDepth);
end

if ~isnumeric(nPops) || nPops<=0
    error('nPops must be a number higher than 0. maxDepth=%d\n',nPops);
end

% search operators
roundTypes = [2*ones(1,nPops) ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 ...
            3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 1 1 1 1 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 ...
            1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];

MAX_SEQ_NUM = 3000;


globalLabelOffset = 0;
for curDepth=0:maxDepth-1
    if REUSE_MODE
        if curDepth<c.depth
            continue;
        end
        
        globalLabelOffset = c.globalLabelOffset;
        
        curPart = c.partArr(c.seqInds,curDepth);
        availCluIDs = unique(curPart(:)');
        
        REUSE_MODE = false;
    elseif curDepth==0
        c = preprocAln(alnMat,maxDepth,1);
        curPart = ones(c.nSeq,1);
        availCluIDs = 1;
        c.heds = heds;
    else
        curPart = c.partArr(c.seqInds,curDepth);
        availCluIDs = unique(curPart(:)');  % note that availCluIDs is in ascending order
    end
    
    fprintf('-------------  curDepth=%d -------------\n\n',curDepth);
    
    if isempty(curPart)
        fprintf('All sequences are clustered. Quit! There are %d layers in total.\n',curDepth);
        break;
    end
    
    curSeqInds = c.seqInds;
    
    localLabelOffset = 0;
    newPart = curPart;

    for i=1:length(availCluIDs)
        cluid = availCluIDs(i);
        
        fprintf('curDepth=%d. clusterID=%d. \n',curDepth,cluid);
        
        if sum(curPart==cluid)<4
            tmpInds = (curPart==cluid);
            c.seqInds = setdiff(c.seqInds,curSeqInds(tmpInds));
            newPart(tmpInds) = globalLabelOffset-1;
            globalLabelOffset = globalLabelOffset - 1;
            continue;
        end
        
        if curDepth==0
            tmpC = c;           
        else
            tmpC = preprocAln(c,curDepth,cluid);            
        end
        
        tmpZ = linkage(tmpC.dist,'complete');
        if tmpC.nSeq>3*nPops
            tmpInitPart = cluster_own(tmpZ,nPops);
        else
            tmpNum = min(round(tmpC.nSeq/2),nPops);
            tmpInitPart = cluster_own(tmpZ,tmpNum);
        end
        
        if tmpC.nSeq<=MAX_SEQ_NUM
            [tmpPart, tmpLogml] = model_search_parallel(tmpC, tmpInitPart, roundTypes);
        else
            tmpMaxDist = max(tmpC.dist);
            tmpMinDist = min(tmpC.dist);
            tmpCutOff = tmpMinDist + 0.05*(tmpMaxDist-tmpMinDist);
            [tmpPart, tmpLogml] = model_search_speed(tmpC, tmpZ, roundTypes, tmpCutOff,length(unique(tmpInitPart)));
        end
        
        tmpUniClusters = unique(tmpPart);
        if length(tmpUniClusters)==1
            % removed the sequences if it can not be splitted
            c.seqInds = setdiff(c.seqInds,tmpC.seqInds);            
            [~,~,tmpInds] = intersect(tmpC.seqInds,curSeqInds);
            newPart(tmpInds) = globalLabelOffset-1;
            globalLabelOffset = globalLabelOffset - 1;
        else
            [~,~,tmpInds] = intersect(tmpC.seqInds,curSeqInds);
            newPart(tmpInds) = tmpPart+localLabelOffset;
            localLabelOffset = localLabelOffset + length(tmpUniClusters);  
        end
    end   % end of inner search in a certain depth
    
    % update ites in c: partArr, data, nSeq (seqInds have been updated already, no need to update dist)
    if curDepth~=0
        c.partArr(:,curDepth+1) = c.partArr(:,curDepth);
        c.depth = c.depth + 1;
    else
        if length(unique(c.partArr(:,1)))>=nPops
            warning('%d clusters detected in the first layer, which is >= nPops=%d. Input a higher nPops.',...
                length(unique(c.partArr(:,1))), nPops);
            partition = c.partArr(:,1);
            logml = tmpLogml;
            return;
        end
        c = rmfield(c,'dist'); % The field 'dist' will not be used anymore
    end
    c.partArr(curSeqInds,curDepth+1) = newPart;
    
    % update items in c: data (seqInds have been updated already)
    [~,~,tmpInds1] = intersect(c.seqInds,curSeqInds);
    tmpData = c.data(tmpInds1,:);
    tmpBases = histc(tmpData,nt2int('ACGT'));
    tmpInds2 = (sum(tmpBases,1)~=max(tmpBases,[],1));
    c.data = tmpData(:,tmpInds2);
    
    % update ites in c: nSeq (seqInds have been updated already)
    c.nSeq = length(c.seqInds);
    
end

c.globalLabelOffset = globalLabelOffset;

% calculate logml and partition
partition = adjustPartArr(c.partArr(:,1:c.depth));
% logml = zeros(1,c.depth);
%for i=1:c.depth
%    logml(i) = checkLogml(c.snpData, partition(:,i));
%end

% dlmwrite('partition.txt',partition,'delimiter','\t');
% dlmwrite('logml.txt',logml,'delimiter','\t','precision', 20);

save(outFile,'c','partition','-v7.3');
dlmwrite(strcat(outFile,'.partition.txt'),partition,'delimiter','\t');
