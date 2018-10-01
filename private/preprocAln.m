function c = preprocAln(c, depth, clusterId)
% c: preprocessed data
% depth: depth of the current hierarchical clustering 
% clusterId: cluster id of the data to be clustered

% Lu Cheng, 18.03.2012

% SNP data stored in 'uint8', c.OrigDist stored in single precision
% calculate pairwise distance within this function
% Lu Cheng, 18.04.2012

% updated 09.08.2013, regarding the SNP site identification and prior specification
% Lu Cheng

if ~isstruct(c)
    alnMat = c;
    clear c
    
    c.nSeq = size(alnMat,1);   % number of sequences that can be further clustered
    c.maxDepth = depth;
    c.partArr = zeros(c.nSeq, c.maxDepth);
    c.partArr(:,1) = 1;
    
    uniqBase = unique(alnMat(:))';
    if ~all(ismember(uniqBase,'-ACGT'))
        error('Unknown base besides ACGT-, unique bases: %s.\n',uniqBase);
    end
    
    alnMat = nt2int(alnMat);
    cntBases = histc(alnMat,nt2int('ACGT'),1);
    c.snpPosition = find(sum(cntBases,1)~=max(cntBases,[],1));               % snp position extracted from the original data
    c.snpData = alnMat(:,c.snpPosition);                                     % snp data extracted from the original data
    
    prior = cntBases(:,c.snpPosition)>0;    % added 09.08.2013, Lu Cheng
    c.prior = prior./(repmat(sum(prior,1),4,1));
    
    clear alnMat cntBases prior
    
%     c.origDist = single(pdist(double(c.snpData),'hamming')*length(c.snpPosition));   % distances between the original sequences
    c.nOrigSeq = c.nSeq;                                                     % number of sequences in the original data
    c.origDist = zeros(1,nchoosek(c.nSeq,2),'single');
    offset = 0;
    for i=1:c.nSeq-1
        tmpN = c.nSeq-i;
        tmpInds = i + (1:tmpN);
        tmpDist = sum(repmat(c.snpData(i,:),tmpN,1)~=c.snpData(tmpInds,:),2);
        c.origDist(offset+(1:tmpN)) = tmpDist';
        offset = offset+tmpN;
    end
    
    clear tmp* offset
        
    c.depth = 1;             % depth of clustering, i.e. current clustering proceeds to c.depth layer
    c.seqInds = 1:c.nSeq;    % sequence indexes that can be further clustered, or available sequences
    c.data = c.snpData;      % snp data correspond to the available sequences
    c.dist = c.origDist;     % distances between the available sequences
    
    return;
end

tmpInds = (c.partArr(c.seqInds,depth)==clusterId);
seqInds = sort(c.seqInds(tmpInds));

if length(seqInds)==1
    c.data = [];
    c.seqInds = seqInds;
    c.dist = [];
    c.nSeq = 1;
    return;
end

data = c.data(tmpInds,:);
cntBases = histc(data,nt2int('ACGT'));
tmpInds = (sum(cntBases,1)~=max(cntBases,[],1));
data = data(:,tmpInds);

prior = cntBases(:,tmpInds)>0;    % added 09.08.2013, Lu Cheng
c.prior = prior./(repmat(sum(prior,1),4,1));

c.data = data;
c.seqInds = seqInds;
c.dist = getDistance(seqInds,c.origDist,c.nOrigSeq);
c.nSeq = length(seqInds);

%----------------------------------------------------------------------------


function [dist2, dind1, dind2] = getDistance(inds2, dist_orig, ninds)
% pick out the distrances between samples in "inds2" from "dist_orig"
% dist_orig specifies the distances of (1,2),(1,3),(1,4)......(ninds-1,ninds)
% Lu Cheng, 22.06.2011

if ~issorted(inds2)
    error('inds2 is not in ascending order!');
end

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
irow = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(irow, 1) = inds2(i);
        apu(irow, 2) = inds2(j);
        irow = irow+1;
    end
end

dind1 = apu(:,1);
dind2 = apu(:,2);

apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = dist_orig(apu);
