function newPlant=slicePlantByZ(plant,zLims,nLay)

    nL=plant.nL;
    cx=plant.cx;
    cy=plant.cy;
    cz=plant.cz;
    parents=plant.parents;
    M=plant.M;
    L=plant.L;
    AX=plant.AX;
    R=plant.R;
    kx=plant.kx;
    kr=plant.kr;
    %also have mSys & nAxes, both are scalar and invariant under slicing

    [inLay,sameLay,crossLay]=linkLayersC(cz,zLims,nL,nLay);
    U=crossUs3(cz,zLims(2:end),inLay,crossLay,nL);
    

    [parents,L,M,nL,cx,cy,cz,R,AX,kr,kx]=breakSegmentsToLayers(U,parents,L,M,nL,cx,cy,cz,R,AX,kr,kx);
    %should import this code here and make/use custom version

    %would have to redo the next one to deal with just coefficients...
    
    [inLay,sameLay,crossLay]=linkLayersC(cz,zLims,nL,nLay); 
    %it may be these are counted as crossing even though they just end on
    %the boundary...
    zEnds=[cz(:,2) sum(cz,2)];

    crossLay=crossLay & ~(ismember(zEnds(:,1),zLims) | ismember(zEnds(:,2),zLims)); %ismember(cz(:,2),zLims);

    if any(crossLay)
        throw(MException('sinkByLay:badLink','Link crosses layer boundary'));
    end
    
    newPlant.nL=nL;
    newPlant.cx=cx;
    newPlant.cy=cy;
    newPlant.cz=cz;
    newPlant.parents=parents;
    newPlant.L=L;
    newPlant.M=M;
    newPlant.mSys=plant.mSys;
    newPlant.nAxes=plant.nAxes;
    newPlant.AX=AX;
    newPlant.R=R;
    newPlant.kr=kr;
    newPlant.kx=kx;
    
    newPlant.inLayer=linkLayersD(zEnds,zLims,nL,nLay);
    
end