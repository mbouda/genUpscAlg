
addpath ~/Documents/restoredDocs/RSAstencil/upscalingJan
%used for checking solutions

myPool=parpool(6); %used in arriving at full solutions

tol=1e-5;  %relative tolerance for soft check
    %for fully determined systems, can do tol=1e-12
    %for overdetermined, get c close to 0, need tol=1e-5
    %regardless, residuals on actual predictions are double-precision roundoffs
resTol=1e-12; %relative tolerance on numerical residuals; 
    %needed for cases where numerical estimation of coefficients rank
    %deficient...
%%
testSet=formTestSet();

nTestSys=size(testSet,2);
collarCond='psiC'; 

for i=1:nTestSys
    testSet(i).nDomLayers=max(testSet(i).inLayer);
    testSet(i).params=testingSetRandParams(testSet(i).parents);
    testSet(i).lyrArch=testSetLyrArch(testSet(i).parents,testSet(i).inLayer,...
        testSet(i).params.nL);
    
    %could these three setup calls be placed outside the loop for speed?
    fldPrp=cell(testSet(i).nDomLayers,1);
    testSet(i).prob=struct('iLinks',fldPrp,'kLayers',fldPrp,'terms',fldPrp,...
        'ints',fldPrp,'bots',fldPrp,'tops',fldPrp,'targ',fldPrp);
    testSet(i).sol=struct('kLayer',fldPrp,'coefs',fldPrp,'vars',fldPrp,'depvar',fldPrp);
    
    %keyboard calls at unresolved points in algorithm can be returned to execution by
    %entering dbcont; should currenly only indicate unfinished parts of code
    
    tic
    for j=1:testSet(i).nDomLayers
        [testSet(i).prob(j),testSet(i).sol(j)]=layerProbSol(j,collarCond,...
            testSet(i).lyrArch,testSet(i).params,testSet(i).parents,testSet(i).inLayer);
    end
    testSet(i).time=toc;
    fullSol=testCaseSol(testSet(i));
    
    testSet(i).check=checkSol(testSet(i),fullSol,tol,resTol);
    testSet(i).passSoft=all(cat(1,testSet(i).check(:).allSoft));
    testSet(i).passRes=all(cat(1,testSet(i).check(:).resTest));
end

result=cat(1,testSet(:).passSoft);
if all(result)
    fprintf(1,'All cases soft-tested successfully to rel. tol=%.3e\n',tol);
else
    fprintf(1,'Failed soft-test in case %d\n',find(~result));
end
result=cat(1,testSet(:).passRes);
if all(result)
    fprintf(1,'All cases tested successfully on residuals to rel. tol=%.3e\n',resTol);
else
    fprintf(1,'Failed residual test in case %d\n',find(~result));
end

times=cat(1,testSet(:).time);

