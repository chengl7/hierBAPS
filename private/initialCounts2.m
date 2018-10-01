function [sumcounts, counts] = initialCounts2(partition, data, npops)
% initialize counts and sumcounts for the initial partition
%    npops: number of populations in the partition

% Lu Cheng, 18.04.2012

% updated 09.08.2013, regarding the SNP site identification and prior specification
% Lu Cheng

global DATA_TYPE;

[~, nLoci] = size(data);

% counts = zeros(5,nLoci,npops);
% sumcounts = zeros(nLoci,npops);

% counts = zeros(5,nLoci,npops,DATA_TYPE);
counts = zeros(4,nLoci,npops,DATA_TYPE); % added 09.08.2013, Lu Cheng
sumcounts = zeros(nLoci,npops,DATA_TYPE);

% letters = nt2int('ACGT-');
letters = nt2int('ACGT'); % added 09.08.2013, Lu Cheng

for i=1:npops
    inds = (partition==i);
    counts(:,:,i) = histc(data(inds,:),letters,1);
    sumcounts(:,i) = sum(counts(:,:,i),1);
end

