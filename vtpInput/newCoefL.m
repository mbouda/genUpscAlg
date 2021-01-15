function [cxNew,cyNew,czNew,l]=newCoefL(plant,u,pcs,iL)

    slopes=[plant.cx(iL,1) plant.cy(iL,1) plant.cz(iL,1)];
    s=sqrt(sum((u*slopes).^2,2));

    l=s(2:end)-s(1:end-1);

    du=u(2:end)-u(1:end-1);
    dx=plant.cx(iL,1)*du;
    dy=plant.cy(iL,1)*du;
    dz=plant.cz(iL,1)*du;

    cxNew=cat(2,dx,u(1:pcs)*plant.cx(iL,1)+plant.cx(iL,2));
    cyNew=cat(2,dy,u(1:pcs)*plant.cy(iL,1)+plant.cy(iL,2));
    czNew=cat(2,dz,u(1:pcs)*plant.cz(iL,1)+plant.cz(iL,2));

end