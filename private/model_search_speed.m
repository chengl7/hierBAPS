function [partition, logml] = model_search_speed(c, Z, roundTypes, cutoff, nMaxPops)
% This function clusters DNA alignment using independent loci model
% c: preprocessed data for the sequence alignment
% partition: initial partition of the individuals
% roundTypes: array of operation types
% cutoff: cutoff to determine sequences belong to a pregroup
% Notice: this function pregroups the squences first, then perform the search
% Lu Cheng
% 11.01.2012

% Data type can be stored in 'uint16' to save memory
% Lu Cheng
% 18.04.2012

% updated 09.08.2013, regarding the SNP site identification and prior specification
% Lu Cheng

interactive = false;

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

nINDS = c.nSeq;
nPOPS = nMaxPops;

% pregroup the sequences, then determine the partition
pgPart = cluster(Z,'cutoff',cutoff,'criterion','distance');
nPregroup = length(unique(pgPart));
if nPregroup<=nMaxPops
    error('#pregroup: %d, nMaxPops: %d. Please lower the cutoff.\n',nPregroup,nMaxPops);
end

pregroups = cell(nPregroup,1);
pgSize = zeros(nPregroup,1);
for i=1:nPregroup
    pregroups{i} = find(pgPart==i);
    pgSize(i) = length(pregroups{i});
end

pgDist = zeros(nchoosek(nPregroup,2),1);
tmp = 1;
for i=1:nPregroup-1
    for j=i+1:nPregroup
        tmpDist = getDistance1(pregroups{i},pregroups{j}, c.dist, nINDS);
        pgDist(tmp)=mean(tmpDist);
        tmp = tmp+1;
    end
end
clear tmp tmpDist

pgZ = linkage(pgDist(:)','complete');
initPart = cluster(pgZ,'maxclust',nPOPS);
partition = zeros(nINDS,1);
for i=1:nPregroup
    partition(pregroups{i}) = initPart(i);
end
clear i pgZ initPart

% Initialize DATA_TYPE,PARTITION, **_COUNTS, SUM_**_COUNTS, alnMat

DATA_TYPE = 'uint16'; % supports up to 65535 input sequences
if c.nSeq>intmax(DATA_TYPE)
    error('DATA_TYPE is %s, maxNum=%d is less than #sequences=%d\n!',...
        DATA_TYPE,intmax(DATA_TYPE),c.nSeq);
end

[sumCounts, counts] = initialCounts2(partition, c.data, nPOPS);

% tmp = sum(counts,3)>0;
% sumTmp = sum(tmp,1);
% tmpPrior = repmat(1./sumTmp,5,1).*tmp;
% priorPar = single(repmat(tmpPrior,[1 1 nPOPS]));

priorPar = single(repmat(c.prior,[1 1 nPOPS])); % 09.08.2013

COUNTS = counts; SUM_COUNTS = sumCounts;

PRIOR_PAR = priorPar;
PARTITION = partition;

clear partition counts sumCounts tmp* sumTmp priorPar;

% Initialize LOGML_TABLE:
LOGML_TABLE = zeros(nPOPS,1);
updateLogmlTable(1:nPOPS);

REMOVAL_DIFFERENCE = zeros(nINDS,1);
REMOVAL_DIFFERENCE(:,:) = nan;
ADDITION_DIFFERENCE = zeros(nINDS,nPOPS);
ADDITION_DIFFERENCE(:,:) = nan;
JOIN_DIFFERENCE = zeros(nPOPS, nPOPS);
JOIN_DIFFERENCE(:,:) = nan;

% ***********Doc:********************
% REMOVAL_DIFFERENCE(ind) tells the change in logml if ind is removed from
% its cluster. nan, if the cluster has changed, since the value was last
% calculated.
% 
% ADDITION_DIFFERENCE(ind, pop) tells the change in logml if ind is added
% to cluster pop. nan, if the cluster has changed since the value was last
% calculated. Always nan, if pop is ind's own cluster.
%
% JOIN_DIFFERENCE(pop1,pop2) = tells the change in logml if pop1 and pop2
% are combined. nan, if either cluster has changed since the value was last
% calculated.
% ***********Doc end*****************

logml = computeTotalLogml();  % uniform K prior

disp('The beginning:');
% disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logml)]);
% disp(' ');

% START SEARCH OF THE BEST PARTITION:

vipu = zeros(1,14); 
if interactive
    roundTypes = input('Input steps: ');
    if ischar(roundTypes), roundTypes = str2num(roundTypes); end
end
ready = 0;


while ready ~= 1

%     disp(['Performing steps: ' num2str(roundTypes)]);

    for n = 1:length(roundTypes)
        round = roundTypes(n);
        moveCounter = 0;

        if  round==1 && vipu(1)==0  % move an individual to another population
            
            pgInds = getMoveInds(pgPart,pgDist,nPregroup); % get pregroup inds to be moved

            for pgind = pgInds(:)'
%                 inds = cell2mat(pregroups(pgInds));
                tmpInds = pregroups{pgind};
                tmpChanges = calcLogmlChanges(tmpInds, c.data, logml);

                [maxChange, maxIndex] = max(tmpChanges);
                if maxChange>1e-5
                    updateGlobalVariables(tmpInds, maxIndex, c.data);
%                     fprintf('moving from %d to %d.\n',PARTITION(ind),maxIndex)
                    logml = computeTotalLogml();
                    moveCounter = moveCounter+length(tmpInds);
                    vipu = zeros(1,14);
                end
            end
            if moveCounter==0, vipu(1)=1; end
            disp(['Step 1: ' num2str(moveCounter) ' individuals were moved.']);
            
        elseif round==2 && vipu(2)==0  % join two populations

            update_join_difference(c.data, logml);
            [maxChange, aux] = max(JOIN_DIFFERENCE(:));
            [i1, i2] = ind2sub([nPOPS,nPOPS],aux);

            if maxChange>1e-5
                tmpInds = find(PARTITION==i1);
                updateGlobalVariables(tmpInds, i2, c.data);
                logml = computeTotalLogml();

                disp(['Step 2: Clusters ' num2str(i1) ' and ' num2str(i2) ' combined.']);
                vipu = zeros(1,14);
            else 
                disp('Step 2: no changes.');
                vipu(2)=1;
            end
        elseif ismember(round, 3:4) && vipu(round)==0  % Split a population, and move one subpopulation to another population

            pops = randperm(nPOPS);

            splitFlags = zeros(nPOPS,1);
            for pop = pops(:)'

                maxChange = 0;
                indsToBeMoved = [];

                inds2 = find(PARTITION==pop);
                pgInds2 = unique(pgPart(inds2));
                nPgInds2 = length(unique(pgPart(inds2)));
                if nPgInds2>4

                    if round==3
                        dist3 = getDistance(pgInds2,pgDist,nPregroup);
                        npops2 = min(20, floor(nPgInds2 / 5));
                    elseif round==4
                        dist3 = getDistance(pgInds2,pgDist,nPregroup);
                        npops2 = 2;
                    end

                    Z3 = linkage(dist3(:)','complete');
                    T3 = cluster(Z3, 'maxclust', npops2);

                    for i = 1:npops2
                        indsX = pgInds2(T3==i);
                        indsX = cell2mat(pregroups(indsX));
                        tmpChanges = calcLogmlChanges(indsX, c.data, logml);
                        [tmpMaxChange, tmpMaxPop] = max(tmpChanges);
                        if tmpMaxChange>maxChange
                            maxChange = tmpMaxChange;
                            % i1 = pop;
                            i2 = tmpMaxPop;
                            indsToBeMoved = indsX;
                        end
                    end
                    if maxChange>1e-5
                        updateGlobalVariables(indsToBeMoved, i2, c.data);
                        logml = computeTotalLogml;
                        splitFlags(pop)=1;
                    end
                end
            end
            if any(splitFlags)
                disp(['Step ' num2str(round) ': ' num2str(sum(splitFlags)) ' populations were split.']);
                vipu = zeros(1,14);
            else
                disp(['Step ' num2str(round) ': no changes.']);
                vipu(round)=1;
            end
        end        
    end

    if interactive
        roundTypes = input('Input extra steps: ');
        if ischar(roundTypes), roundTypes = str2num(roundTypes); end
    else
        roundTypes = [];
    end

    if isempty(roundTypes)
        ready = 1;
    end
end

nPOPS = rmEmptyPopulation();
partition = PARTITION;

% disp(' ');
disp('BEST PARTITION: ');
% disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logml)]);
disp(' ');


%----------------------------------------------------------------------------


function [dist2, dind1, dind2] = getDistance(inds2, origDist, ninds)
% pick out the distrances between samples in "inds2" from "origDist"
% origDist specifies the distances of (1,2),(1,3),(1,4)......(ninds-1,ninds)
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
dist2 = origDist(apu);

function dist = getDistance1(inds1,inds2, origDist, ninds)
% pick out the distrances between samples in "inds1" and "inds2" from "origDist"
% "inds1" and "ind2" should have no common elements
% origDist specifies the distances of (1,2),(1,3),(1,4)......(ninds-1,ninds)
% Lu Cheng, 22.06.2011

if any(ismember(inds1,inds2))
    error('inds1 and inds2 have common elements.');
end

ninds1 = length(inds1);
ninds2 = length(inds2);
apu = zeros(ninds1*ninds2,2);

[x1, x2] = meshgrid(inds1,inds2);
apu(:,1) = min(x1(:),x2(:));
apu(:,2) = max(x1(:),x2(:));

apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist = origDist(apu);

%---------------------------------------------------------------


function inds = getMoveInds(pgPart, pgDist, nPregroup)
% get pregroup indexs to be moved to another cluster
% we always take the 35% pregroups of each cluster which are most distant
% to each other
% Lu Cheng, 22.06.2011

global PARTITION;

pops = unique(PARTITION);
inds = [];

for tmpPop = pops(:)'
    tmpInds = unique(pgPart(PARTITION==tmpPop));

    if(length(tmpInds)<20)
        inds = [inds tmpInds(:)']; %#ok<AGROW>
        continue;
    end
    
    [tmpDist, dind1, dind2] = getDistance(tmpInds,pgDist,nPregroup);
    tmpVal = quantile(tmpDist,0.65);
    tmpInds2 = find(tmpDist>tmpVal);
    tmpInds3 = union(unique(dind1(tmpInds2)), unique(dind2(tmpInds2)));
    inds = [inds tmpInds3(:)']; %#ok<AGROW>
end










