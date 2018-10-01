function exData(inFile,type)
% read the input sequence file, and transform it to alignment matrix
% Lu Cheng
% 12.09.2012

% modified on 16.1.2015, added the '-v7.3' flag to save commands
% Lu Cheng

if strcmp(inFile(end-3:end),'.xls') || strcmp(type,'xls')
    [~,~,rawdata] = xlsread(inFile);
    heds = cellfun(@num2str,rawdata(2:end,1),'UniformOutput',false);
    data = rawdata(2:end,2:end);
    alnMat = cell2mat(data);
    alnMat = upper(alnMat);
save('seqs.mat','alnMat','heds','-v7.3');
    return;

elseif strcmp(inFile(end-5:end),'.fasta') || strcmp(type,'fasta')
    [heds, seqs]=readFasta(inFile);
    for i=1:length(seqs)
        inds = ~ismember(seqs{i},'ACGTacgt');
        if sum(inds)>0
            seqs{i}(inds)='-';
        end
    end
    alnMat = cell2mat(seqs);
    alnMat = upper(alnMat);
save('seqs.mat','alnMat','heds','-v7.3');
    return;
else
    fprintf('Unknown input file type! Should be .xls or .fasta. \n');
    return;
end
