function [plant,zMin,zLims,dz]=importPlant(fileName,nLayInit,kr,kx,b)

    vtpRoot=readRootVTP(fileName);
    vtPlant=formPlantFromVTP(vtpRoot);
    
    zMin=floor(min(sum(vtPlant.cz,2)))/100;
    dz=zMin/nLayInit;
    zLims=(0:dz:zMin)';
    
    vtPlant=adjustVTPlant2(vtPlant);
    
    [vtPlant.kr,vtPlant.kx]=genHydroConst(kr,kx,b,vtPlant.R);
    %should probably re-categorise  some of those fields into subfields of
    %group-heading fields (geom...)
    
    posStart=vtPlant.cz(:,2)>0;
    posEnd=sum(vtPlant.cz,2)>0;
    anyPos=posStart | posEnd;
    
    if any(anyPos)   %if any segments stick out of the ground
        maxZ=max(cat(1,vtPlant.cz(:,2),sum(vtPlant.cz,2)));
        dz=(zMin-maxZ)/nLayInit;
        zLims=(maxZ:dz:zMin)'; %add maximum Z vertex to the layer limits
        vtPlant.kr(anyPos)=0; %set kr any protruding segments to 0
    end            
    
    slPlant=slicePlantByZ(vtPlant,zLims,nLayInit);
    plant=forceTriJuncs(slPlant); 
    
    plant.b=repmat(b,[plant.nL 1]);
    plant.params=formPlantParams(plant);

end