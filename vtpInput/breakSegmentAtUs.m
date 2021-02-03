function newPlant=breakSegmentAtUs(nuBrk,iL,plant,badPar)

%works only for linear case
     newPlant=plant;
    
    %actually, need to find set of 'descendant' links,
    %then can nicely add +1 every time we split off a new link 
    
%    desc=allDescendents(iL,plant.parents);
    
    u=cat(1,0,nuBrk,1);
    pcs=size(u,1)-1;
    
    if iL==plant.parents(iL)+1
        newPar=plant.parents(iL)+(0:pcs-1)';
    else
        newPar=cat(1,plant.parents(iL),iL+(0:pcs-2)');
    end
    
    [cxNew,cyNew,czNew,l]=newCoefL(plant,u,pcs,iL);
    
    %now, need to inject new values into arrays instead of single old
    %values
    
    bef=(1:iL-1)';
    aft=(iL+1:plant.nL)';
    
    newPlant.cx=cat(1,newPlant.cx(bef,:),cxNew,newPlant.cx(aft,:));
    newPlant.cy=cat(1,newPlant.cy(bef,:),cyNew,newPlant.cy(aft,:));
    newPlant.cz=cat(1,newPlant.cz(bef,:),czNew,newPlant.cz(aft,:));
    newPlant.L=cat(1,newPlant.L(bef),l,newPlant.L(aft));
    newPlant.M=cat(1,newPlant.M(bef),repmat(newPlant.M(iL),[pcs,1]),newPlant.M(aft));
    
    newPlant.AX=cat(1,newPlant.AX(bef),repmat(newPlant.AX(iL),[pcs,1]),newPlant.AX(aft));
    newPlant.R=cat(1,newPlant.R(bef),repmat(newPlant.R(iL),[pcs,1]),newPlant.R(aft));
    newPlant.kr=cat(1,newPlant.kr(bef),repmat(newPlant.kr(iL),[pcs,1]),newPlant.kr(aft));
    newPlant.kx=cat(1,newPlant.kx(bef),repmat(newPlant.kx(iL),[pcs,1]),newPlant.kx(aft));
    newPlant.inLayer=cat(1,newPlant.inLayer(bef),repmat(newPlant.inLayer(iL),[pcs,1]),newPlant.inLayer(aft));
    
    %problem here is if badpar==0& iL=1
    if badPar==0
        newPlant.parents(plant.parents>=iL)=plant.parents(plant.parents>=iL)+pcs-1;
    else
        newPlant.parents(plant.parents>iL)=plant.parents(plant.parents>iL)+pcs-1;
    end
    
    newPlant.parents=cat(1,newPlant.parents(bef),newPar,newPlant.parents(aft));
    mainD=find(plant.parents==badPar,1,'first')+pcs-1;
    newPlant.parents(mainD)=badPar+pcs-1; %this should add the 'main daughter' to the last of the new pieces

    newPlant.nL=size(newPlant.M,1);

end