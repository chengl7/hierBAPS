function [partition, logml] = model_search_parallel(c, partition, roundTypes)
% This function clusters DNA alignment using independent loci model
% c: preprocessed data for the sequence alignment
% partition: initial partition of the individuals
% roundTypes: array of operation types

% Lu Cheng
% 07.02.2012

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

nPOPS = length(unique(partition));

% Initialize DATA_TYPE, PARTITION, COUNTS, SUM_COUNTS, PRIOR_PAR, PARTITION
DATA_TYPE = 'uint16'; % supports up to 65535 input sequences
                      % if this value changes, please update also model_search_speed.m, checkLogml.m
if c.nSeq>intmax(DATA_TYPE)
    error('DATA_TYPE is %s, maxNum=%d is less than #sequences=%d\n!',...
        DATA_TYPE,intmax(DATA_TYPE),c.nSeq);
end

[sumCounts, counts] = initialCounts2(partition, c.data, nPOPS);

% tmp = sum(counts,3)>0;
% sumTmp = sum(tmp,1);
% tmpPrior = repmat(1./sumTmp,5,1).*tmp;
% priorPar = single(repmat(tmpPrior,[1 1 nPOPS]));

priorPar = single(repmat(c.prior,[1 1 nPOPS])); % 09.08.2013 Lu Cheng

COUNTS = counts; SUM_COUNTS = sumCounts;

PRIOR_PAR = priorPar;
PARTITION = partition;

clear partition counts sumCounts tmp* sumTmp priorPar;

% Initialize LOGML_TABLE:
nINDS = c.nSeq;
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
            
%             inds = randperm(nINDS); 
            inds = getMoveInds(c.dist,nINDS); % get inds to be moved

            for ind = inds(:)'
                update_difference_tables(ind, c.data, logml);
                tmpDiff = REMOVAL_DIFFERENCE(ind) + ADDITION_DIFFERENCE(ind,:);
                tmpDiff(PARTITION(ind)) = 0;
                [maxChange, maxIndex] = max(tmpDiff);
                if maxChange>1e-5
                    updateGlobalVariables(ind, maxIndex, c.data);
%                     fprintf('moving from %d to %d.\n',PARTITION(ind),maxIndex)
                    logml = computeTotalLogml();
                    moveCounter = moveCounter+1;
                    vipu = zeros(1,14);
                end
            end
            if moveCounter==0, vipu(1)=1; end
            %disp(['Step 1: ' num2str(moveCounter) ' individuals were moved.']);

        elseif round==2 && vipu(2)==0  % join two populations

            update_join_difference(c.data, logml);
            [maxChange, aux] = max(JOIN_DIFFERENCE(:));
            [i1, i2] = ind2sub([nPOPS,nPOPS],aux);

            if maxChange>1e-5
                tmpInds = find(PARTITION==i1);
                updateGlobalVariables(tmpInds, i2, c.data);
                logml = computeTotalLogml();

                %disp(['Step 2: Clusters ' num2str(i1) ' and ' num2str(i2) ' combined.']);
                vipu = zeros(1,14);
            else 
                %disp('Step 2: no changes.');
                vipu(2)=1;
            end
        elseif ismember(round, 3:4) && vipu(round)==0  % Split a population, and move one subpopulation to another population

            pops = randperm(nPOPS);

            splitFlags = zeros(nPOPS,1);
            for pop = pops(:)'

                maxChange = 0;
                indsToBeMoved = [];

                inds2 = find(PARTITION==pop);
                ninds2 = length(inds2);
                if ninds2>4

                    if round==3
                        dist3 = getDistance(inds2, c.dist, nINDS);
                        npops2 = min(20, floor(ninds2 / 5));  %Moneenko osaan jaetaan    
                    elseif round==4
                        dist3 = getDistance(inds2, c.dist, nINDS);
                        npops2 = 2;
                    end

                    Z3 = linkage(dist3);
                    T3 = cluster_own(Z3, npops2);

                    for i = 1:npops2
                        indsX = inds2(T3==i); indsX = indsX';
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
                        logml = computeTotalLogml();
                        splitFlags(pop)=1;
                    end
                end
            end
            if any(splitFlags)
                %disp(['Step ' num2str(round) ': ' num2str(sum(splitFlags)) ' populations were split.']);
                vipu = zeros(1,14);
            else
                %disp(['Step ' num2str(round) ': no changes.']);
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
%disp(['Partition: ' num2str(PARTITION')]);
disp(['Nclusters: ' num2str(length(unique(PARTITION)))]);
disp(['Log(ml*prior): ' num2str(logml)]);
disp(' ');

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


%---------------------------------------------------------------


function inds = getMoveInds(dist_orig, ninds)
% get individual indexs to be moved to another cluster
% we always take the 30% individuals of each cluster which are most distant
% to each other
% Lu Cheng, 25.05.2011

global PARTITION;

pops = unique(PARTITION);
inds = [];

for tmpPop = pops(:)'
    tmpInds = find(PARTITION == tmpPop)';
       
    if(length(tmpInds)<20)
        inds = [inds tmpInds(:)']; %#ok<AGROW>
        continue;
    end
    
    [tmpDist, dind1, dind2] = getDistance(tmpInds,dist_orig,ninds);
    tmpSDist = sort(tmpDist,'Descend');
    tmpInds2 = find(tmpDist>tmpSDist(round(length(tmpSDist)*0.3)));
    tmpInds3 = union(unique(dind1(tmpInds2)), unique(dind2(tmpInds2)));
    inds = [inds tmpInds3(:)']; %#ok<AGROW>
end










