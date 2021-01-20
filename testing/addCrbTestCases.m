function testSet=addCrbTestCases(testSet,dirName,fileName)

    dat=load(strcat(dirName,fileName));
    testSet=cat(2,testSet,dat.testSet);
end