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
dataDir='./testing/crbTestSet/';

nDay=13; 
rsaFile=sprintf('RLab_210119_Pisum_sativum_a_Pagès_2014_%ddenni_simulace.vtp',nDay);

%rsaFile='RLab_210115_workshop.vtp';

kx=5e-5;
kr=1.5e-10;
b=100e-6;
nLayInit=8;

[plant,zMin,zLims,dz]=importPlant(strcat(dataDir,rsaFile),nLayInit,kr,kx,b); 

collarCond='psiC'; 

%% Upscaling

    plant.nDomLayers=nLayInit;
    plant.lyrArch=setLyrArch(plant.parents,plant.inLayer,plant.nL);
    fldPrp=cell(plant.nDomLayers,1);
    
    plant.prob=struct('iLinks',fldPrp,'kLayers',fldPrp,'terms',fldPrp,...
        'ints',fldPrp,'bots',fldPrp,'tops',fldPrp,'targ',fldPrp);
    plant.sol=struct('kLayer',fldPrp,'coefs',fldPrp,'vars',fldPrp,'depvar',fldPrp);
    
    %in Pisum sativum 13-day, j=3: hangning sector in the bottom subdomain
    %layer,  comes from below, ultimately suspended from targ;
    %current way of searching for how to attach it leaves targ as parent
    %and looks for landing on sibling, failing to attach
    
    tic
    for j=1:plant.nDomLayers
        [plant.prob(j),plant.sol(j)]=layerProbSol(j,collarCond,...
            plant.lyrArch,plant.params,plant.parents,plant.inLayer);
    end
    plant.time=toc;
    
  
    
    plant.params.b=plant.b;
    plant.params.r=plant.R;
    plant.params.kr=plant.kr;
    plant.params.L=plant.L;
    tic
    testSol=fullSol(plant);
    toc
    
    plant.check=checkSol(plant,testSol,tol,resTol);
    plant.passSoft=all(cat(1,plant.check(:).allSoft));
    plant.passRes=all(cat(1,plant.check(:).resTest));


%% Evaluation
result=cat(1,plant.passSoft);
if all(result)
    fprintf(1,'Plant soft-tested successfully to rel. tol=%.3e\n',tol);
else
    fprintf(1,'Plant failed soft-test\n');
end
result=cat(1,plant.passRes);
if all(result)
    fprintf(1,'Plant tested successfully on residuals to rel. tol=%.3e\n',resTol);
else
    fprintf(1,'Plant failed residual test\n');
end

resMax=cat(1,plant.check(:).maxRes);

times=cat(1,plant.time); %times of upscaling, not solution; those are trivial.


%% Outputs

% Should probably come out in some recognisable form, e.g. C matrix from Jan's
% approach... in this case, should look how that was constructed for
% simple upscaled model in the case presented at LLN