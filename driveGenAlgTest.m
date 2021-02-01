
addpath ./algo
addpath ./testing
addpath ./fullSol
addpath ./vtpInput
myPool=parpool(6); %used in arriving at full solutions

tol=1e-5;  %relative tolerance for soft check
    %for fully determined systems, can do tol=1e-12
    %for overdetermined, get c close to 0, need tol=1e-5
    %regardless, residuals on actual predictions are double-precision roundoffs
    %system 9 fails softtest because of the overdetermined numerical
    %solution shape, not because of the upscaling algorithm; residuals OK
resTol=1e-12; %relative tolerance on numerical residuals; 
    %needed for cases where numerical estimation of coefficients rank
    %deficient...
    %have now registered some cases of residuals ~3e-12, so fails but fine,
    %soft-test is passed in these cases; all pass if re-run: due to random
    %variation in parameters.
%%
testSet=formTestSet();
nSimpleSys=size(testSet,2);
[testSet,crbParams]=addCrbTestCases(testSet,'./testing/crbTestSet/','testSet.mat');

nTestSys=size(testSet,2);
collarCond='psiC'; 

% max residuals elevated but perhaps acceptable; 
%may be improved by pivoting in final system reduction (?)
% i=19, j=4 cca 4e-12 (likely not random, though)
% i=23, j=4 cca 4e-9
% i=27, j=2 cca 2e-10 (13-day P. sativum) 

%Currently choose just one extra layerEq to use;
    %but don't have a good rule for it... may need to update
%have a way to add equations just before pivoting, but again, may need improvement.

%present last test system (Pisum sativum, 16-day old) is not practical for tests:
%takes way too long and does not appear to add any issues

%should perhaps upgrade use of parallelism in testing, for more elaborate
%cases (thousands of network segments)

for i=1:nTestSsys-1
    testSet(i).nDomLayers=max(testSet(i).inLayer);
    if i<=nSimpleSys
        testSet(i).params=testingSetRandParams(testSet(i).parents);  
    else
        %testSet(i).params=addCrbParams(testSet(i),crbParams(i-nSimpleSys));
        testSet(i).params=crbParams(i-nSimpleSys);
    end
    
    testSet(i).lyrArch=setLyrArch(testSet(i).parents,testSet(i).inLayer,...
        testSet(i).params.nL);
    
    %could these three setup calls be placed outside the loop for speed?
    fldPrp=cell(testSet(i).nDomLayers,1);
    testSet(i).prob=struct('iLinks',fldPrp,'kLayers',fldPrp,'terms',fldPrp,...
        'ints',fldPrp,'bots',fldPrp,'tops',fldPrp,'targ',fldPrp);
    testSet(i).sol=struct('kLayer',fldPrp,'coefs',fldPrp,'vars',fldPrp,'depvar',fldPrp);
    
    %keyboard calls should currenly only indicate unfinished parts of code
    %noted in issues; can be returned to execution by entering dbcont.
    
    tic
    for j=1:testSet(i).nDomLayers
        [testSet(i).prob(j),testSet(i).sol(j)]=layerProbSol(j,collarCond,...
            testSet(i).lyrArch,testSet(i).params,testSet(i).parents,testSet(i).inLayer);
    end
    testSet(i).time=toc;
    testSol=fullSol(testSet(i));
    
    testSet(i).check=checkSol(testSet(i),testSol,tol,resTol);
    testSet(i).passSoft=all(cat(1,testSet(i).check(:).allSoft));
    testSet(i).passRes=all(cat(1,testSet(i).check(:).resTest));
end

%% Evaluation
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

maxRes=cat(1,testSet(17).check(:).maxRes);

%looks like currently some residuals grow out of hand in the crb cases
    %appears to be due to random assignment of parameters -- with params
    %from crb functions, error <1e-12

%seleciton of solution/elimination variables at end to be addressed...
    %tries to minimise the number of coefficients, though not sure if well...
        %if had more of them, they might be smaller, reducing error?
        
times=cat(1,testSet(:).time); %times of upscaling, not solution; those are trivial.

