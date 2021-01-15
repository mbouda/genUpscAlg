function [parents,L,M,nL,cx,cy,cz,r,AX,kr,kx]=breakSegmentsToLayers(nuBrk,parents,L,M,nL,cx,cy,cz,r,AX,kr,kx)

%works only for linear case
 
    j=1;
    parBrk=[];
    LBrk=[];
    MBrk=[];
    cxBrk=[];
    cyBrk=[];
    czBrk=[];
    rBrk=[];
    axBrk=[];
    kxBrk=[];
    krBrk=[];
    
    midWf=zeros(nL,1);
    for i=1:nL
        if any(nuBrk{i}~=0 & nuBrk{i}~=1) %this misses 0s (correctly) but collects 1s ...
            u=unique(cat(1,0,nuBrk{i},1));
            pcs=size(u,1)-1;
            
            slopes=[cx(i,1) cy(i,1) cz(i,1)];           
            s=sqrt(sum((u*slopes).^2,2));
            
            l=s(2:end)-s(1:end-1);
            
            du=u(2:end)-u(1:end-1);
            dx=cx(i,1)*du;
            dy=cy(i,1)*du;
            dz=cz(i,1)*du;
            
            cxNew=cat(2,dx,u(1:pcs)*cx(i,1)+cx(i,2));
            cyNew=cat(2,dy,u(1:pcs)*cy(i,1)+cy(i,2));
            czNew=cat(2,dz,u(1:pcs)*cz(i,1)+cz(i,2));
            
            
            parBrk(j)=midWf(i);
            LBrk(j)=l(1);
            MBrk(j)=M(i);
            cxBrk(j,:)=cxNew(1,:);
            cyBrk(j,:)=cyNew(1,:);
            czBrk(j,:)=czNew(1,:);
            rBrk(j)=r(i);
            axBrk(j)=AX(i);
            krBrk(j)=kr(i);
            kxBrk(j)=kx(i);
            
            j=j+1;
            for k=2:pcs
                parBrk(j)=j-1;
                LBrk(j)=l(k);
                MBrk(j)=M(i);
                cxBrk(j,:)=cxNew(k,:);
                cyBrk(j,:)=cyNew(k,:);
                czBrk(j,:)=czNew(k,:);
                rBrk(j)=r(i);
                axBrk(j)=AX(i);
                krBrk(j)=kr(i);
                kxBrk(j)=kx(i);
                j=j+1;
            end
            midWf(parents==i)=j-1;
        else
            parBrk(j)=midWf(i);
            LBrk(j)=L(i);
            MBrk(j)=M(i);
            cxBrk(j,:)=cx(i,:);
            cyBrk(j,:)=cy(i,:);
            czBrk(j,:)=cz(i,:);
            rBrk(j)=r(i);
            axBrk(j)=AX(i);
            krBrk(j)=kr(i);
            kxBrk(j)=kx(i);
            j=j+1;
            midWf(parents==i)=j-1;
        end
        
    end
    
    parents=parBrk';
    L=LBrk';
    M=MBrk';
    r=rBrk';
    cx=cxBrk;
    cy=cyBrk;
    cz=czBrk;
    AX=axBrk';
    kr=krBrk';
    kx=kxBrk';
    
    nL=size(M,1);
    


end