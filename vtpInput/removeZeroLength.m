function plant=removeZeroLength(plant)

    isZero=plant.L==0;
    J=find(isZero)';
    for j=J
        dtrs=plant.parents==j;
        par=plant.parents(j);
        plant.parents(dtrs)=par;
        plant.parents(plant.parents>j)=plant.parents(plant.parents>j)-1;
    end
    
    plant.cx(J,:)=[];
    plant.cy(J,:)=[];
    plant.cz(J,:)=[];
    
    plant.parents(J)=[];
    plant.L(J)=[];
    plant.M(J)=[];
    plant.AX(J)=[];
    plant.R(J)=[];
    
    plant.nL=plant.nL-sum(isZero);
end