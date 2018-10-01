function [counts sumcounts] = computeDiffInCounts(rows, data)
% calculate the counts of the given rows of the data (ninds*nLoci)

% Lu Cheng, 18.04.2012

% updated 09.08.2013, regarding the SNP site identification and prior specification
% Lu Cheng

global DATA_TYPE;

letters = nt2int('ACGT');
counts = cast(histc(data(rows,:),letters,1),DATA_TYPE);
sumcounts = cast(sum(counts,1)',DATA_TYPE);
