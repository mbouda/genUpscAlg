function [testSet,crbParams]=addCrbTestCases(testSet,dirName,fileName)


    dat=load(strcat(dirName,fileName));
    crbParams=cat(1,dat.testSet(:).params);
    
    dat.testSet=rmfield(dat.testSet,'params');
    testSet=cat(2,testSet,dat.testSet);
    
    
end