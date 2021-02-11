function [zMin,zLims,dz]=setzLims(plant,nLayInit,z0,k)

    lowestZ=min(cat(1,plant.cz(:,end),sum(plant.cz,2)));
    zLims=-inf(2,1);
    c=0.1;
    while sum(k*zLims<lowestZ)>1
        c=10*c;
        
        zMin=floor(c*lowestZ)/(c*k);

        dz=zMin/nLayInit;
        zLims=(z0:dz:zMin)';
    end
end