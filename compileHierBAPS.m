fileDir = 'files';
outDir = 'hierBAPS_package';

if ~exist(outDir,'dir')
    mkdir(outDir);
else
    rmdir(outDir,'s');
    mkdir(outDir);
end

mcc('-m','exData.m','-d',outDir);
mcc('-m','hierBAPS.m','-d',outDir);
mcc('-m','drawSnpMat.m','-d',outDir);

str =  computer;
if strcmp(str,'GLNXA64')
    delete([outDir filesep '*.c']);
    delete([outDir filesep '*.prj']);
    delete([outDir filesep '*.log']);
    delete([outDir filesep '*.sh']);

    movefile([outDir filesep 'readme.txt'],[outDir filesep 'mcr_readme.txt']);
    
    sysName = 'linux';
    
    copyfile([fileDir filesep 'hierBAPS_' sysName '.sh'], [outDir filesep 'hierBAPS.sh']);
    copyfile([fileDir filesep 'readme_' sysName '.txt'], [outDir filesep]);
    copyfile([fileDir filesep 'seqs.fa'], [outDir filesep]);
        
    zip(['hierBAPS_' sysName '_64bit.zip'],outDir);
    
elseif strcmp(str,'MACI64')
    delete([outDir filesep '*.log']);
    delete([outDir filesep '*.sh']);
    
    movefile([outDir filesep 'readme.txt'],[outDir filesep 'mcr_readme.txt']);
    sysName = 'mac';
    
    copyfile([fileDir filesep 'hierBAPS_' sysName '.sh'], [outDir filesep 'hierBAPS.sh']);
    copyfile([fileDir filesep 'readme_' sysName '.txt'], [outDir filesep]);
    copyfile([fileDir filesep 'seqs.fa'], [outDir filesep]);

    zip(['hierBAPS_' sysName '_64bit.zip'],outDir);

elseif strcmp(str,'PCWIN64')
    delete([outDir filesep '*.log']);
    movefile([outDir filesep 'readme.txt'],[outDir filesep 'mcr_readme.txt']);
    
    sysName = 'win';
    
    copyfile([fileDir filesep 'readme_' sysName '.txt'], [outDir filesep]);
    copyfile([fileDir filesep 'seqs.fa'], [outDir filesep]);
    copyfile([fileDir filesep 'seqs.xls'], [outDir filesep]);
        
    zip(['hierBAPS_' sysName '_64bit.zip'],outDir);
else
    error('%s is not supported yet.\n',str);
end