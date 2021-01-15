% This script is used to drive the upscaling process on an RSA created by
% CPlantBox, stored in a *.vtp file.

%% Preliminaries
addpath ./algo ./fullSol ./testing ./vtpInput

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
%% File input, preprocessing
dataDir='/run/media/mbouda/OS/Users/Martin/Documents/graphics/';
rsaFile='Triticum_aestivum_LAB201006.vtp';

kx=5e-5;
kr=1.5e-10;
b=100e-6;
nLayInit=8;

[plant,zMin,zLims,dz]=importPlant(strcat(dataDir,rsaFile),nLayInit,kr,kx,b); %uses same system as in "upscalePlant" directory

collarCond='psiC'; 

%% Upscaling


    plant.nDomLayers=nLayInit;
    plant.lyrArch=setLyrArch(plant.parents,plant.inLayer,plant.nL);
    fldPrp=cell(plant.nDomLayers,1);
    plant.prob=struct('iLinks',fldPrp,'kLayers',fldPrp,'terms',fldPrp,...
        'ints',fldPrp,'bots',fldPrp,'tops',fldPrp,'targ',fldPrp);
    plant.sol=struct('kLayer',fldPrp,'coefs',fldPrp,'vars',fldPrp,'depvar',fldPrp);
    
    tic
    for j=1:plant.nDomLayers
        [plant.prob(j),plant.sol(j)]=layerProbSol(j,collarCond,...
            plant.lyrArch,plant.params,plant.parents,plant.inLayer);
    end
    plant.time=toc;
    testSol=fullSol(plant);
    
    plant.check=checkSol(plant,testSol,tol,resTol);
    plant.passSoft=all(cat(1,plant.check(:).allSoft));
    plant.passRes=all(cat(1,plant.check(:).resTest));


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

times=cat(1,testSet(:).time); %times of upscaling, not solution; those are trivial.


%% Outputs

% Should probably come out in some recognisable form, e.g. C matrix from Jan's
% approach... in this case, should look how that was constructed for
% simple upscaled model in the case presented at LLN