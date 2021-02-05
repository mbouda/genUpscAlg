function plant=removeZeroLength(plant)

    isZero=plant.L==0;
    J=fliplr(find(isZero)');
    for j=J
        dtrs=plant.parents==j;
        par=plant.parents(j);
        plant.parents(dtrs)=par;
        plant.parents(plant.parents>j)=plant.parents(plant.parents>j)-1;
        
        plant.cx(j,:)=[];
        plant.cy(j,:)=[];
        plant.cz(j,:)=[];
    
        plant.parents(j)=[];
        plant.L(j)=[];
        plant.M(j)=[];
        plant.AX(j)=[];
        plant.R(j)=[];        
    end
    
    plant.nL=plant.nL-sum(isZero);
end