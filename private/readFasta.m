function [heds seqs]=readFasta(inFile)
% This function reads large fasta files
% Lu Cheng
% 12.09.2012

nSeq=0;
nLine=0;
lens=[];
seqLineNum = [];

fid = fopen(inFile);
tLine = strtrim(fgetl(fid));

while ischar(tLine)
    tLine = strtrim(tLine);
    nLine=nLine+1;
    if isempty(tLine)
        seqLineNum = seqLineNum + 1;
    elseif tLine(1)=='>'
        lens = [lens seqLineNum];
        nSeq=nSeq+1;
        seqLineNum = 0;
    else
        seqLineNum = seqLineNum + 1;
    end
    tLine = fgetl(fid);
end
lens = [lens seqLineNum];
fclose(fid);

heds = cell(nSeq,1);
seqs = cell(nSeq,1);

fid = fopen(inFile);
for i=1:nSeq
    % read header
    tmpLine = fgetl(fid);
    heds{i}=tmpLine(2:end);
    
    % read sequence
    tmpLines = cell(1,lens(i));
    for j=1:lens(i)
        tmpLines{j}=fgetl(fid);
    end
    seqs{i}=cell2mat(tmpLines);
end
