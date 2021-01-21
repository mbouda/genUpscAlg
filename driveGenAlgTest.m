
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
testSet=addCrbTestCases(testSet,'./testing/crbTestSet/','testSet.mat');

nTestSys=size(testSet,2);
collarCond='psiC'; 

%newest changes appear to have stranded link 7 in case 1, layer 5
    %redefine problem or hooking?
    %or is it that we should close target (8)... that way, can use 2 Eqs to
    %elim the 2 remaining vars
        %looks like getPars does not identify link 8 as par of link 9 ...
        %not sure why?

for i=1:13 %nTestSys
    testSet(i).nDomLayers=max(testSet(i).inLayer);
    testSet(i).params=testingSetRandParams(testSet(i).parents);
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

%looks like currently some residuals on order 1e-10 in the crb cases,
%likely due to seleciton of solution/elimination variables at end; to be
%addressed...

times=cat(1,testSet(:).time); %times of upscaling, not solution; those are trivial.

