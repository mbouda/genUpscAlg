% This script is used to upscaling and check solutions for CPlantBox *.vtp
% files in a given directory

%% Preliminaries
addpath ./algo ./fullSol ./testing ./vtpInput

tol=1e-5;  
resTol=1e-12; 
%% File input, preprocessing
dataDir='./testing/crbTestSet/RLab_210209_plantSet/';

files=dir(dataDir);
fileNames=cat(2,{files(:).name})';
fileNames=fileNames(startsWith(fileNames,'RLab'));
nFiles=size(fileNames,1);

kx=5e-5;
kr=1.5e-10;
b=100e-6;
nLayInit=8;
collarCond='psiC'; 

myPool=parpool(8); 

result=zeros(nLayInit,nFiles);

for i=1:nFiles

    fprintf(1,'Starting %s \n',fileNames{i});
    
    try
        [plant,zMin,zLims,dz]=importPlant(strcat(dataDir,fileNames{i}),nLayInit,kr,kx,b); 
    catch
        warning('Failed to import, skipping...\n','noImp')
        result(:,i)=-1;
        continue
    end

    if any(plant.params.c2==0)
        warning('c2==0 for some links, leading to NaN solns\n likely due to L<1e-8','badL')
    end

    plant.nDomLayers=nLayInit;
    lyrArch=setLyrArch(plant.parents,plant.inLayer,plant.nL);

    params=plant.params;
    parents=plant.parents;
    inLayer=plant.inLayer;
    
    parfor j=1:plant.nDomLayers
        try
            [plant(j).prob,plant(j).sol]=layerProbSol(j,collarCond,...
                lyrArch,params,parents,inLayer);
        catch
            warning(sprintf('Failed to uscale layer %d\n',j),'noUpsc');
            result(j,i)=-2;
            
            fldPrp=cell(1,1);
            plant(j).prob=struct('iLinks',fldPrp,'kLayers',fldPrp,'terms',fldPrp,...
                'ints',fldPrp,'bots',fldPrp,'tops',fldPrp,'targ',fldPrp,'eqs',fldPrp);
            plant(j).sol=struct('kLayer',fldPrp,'coefs',fldPrp,'vars',fldPrp,'depvar',fldPrp);
    
            plant(j).sol.depvar=sprintf('psiXBar%d',j);
            js=setdiff([1 4 8],j);
            plant(j).sol.vars={sprintf('psiL%d',j), sprintf('psiL%d',js(1)),...
                sprintf('psiXBar%d',js(1)), sprintf('psiXBar%d',js(2))};
            plant(j).sol.coefs=zeros(1,4);
        end
    end
    
    for j=2:plant(1).nDomLayers
        plant(1).prob(j)=plant(j).prob;
        plant(1).sol(j)=plant(j).sol;
    end
    for j=plant(1).nDomLayers:-1:2
        plant(j)=[];
    end

    plant.params.b=plant.b;
    plant.params.r=plant.R;
    plant.params.kr=plant.kr;
    plant.params.L=plant.L;

    try
        testSol=fullSol(plant);
    catch
        warning('Failed to find full solution \n','noSol')
        result(:,i)=-3;
        continue
    end

    plant.check=checkSol(plant,testSol,tol,resTol);
    resMax=cat(1,plant.check(:).maxRes);

    result(result(:,i)==0,i)=resMax(result(:,i)==0);
    
end

save('testSweep','fileNames','result');

