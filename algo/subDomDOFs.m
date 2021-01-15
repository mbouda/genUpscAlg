function nDOF=subDomDOFs(top,bot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer)

    %at top:
    xUp=inLayer==top & (nxtLr==(top-1) | prvLr==(top-1));
    nXUp=sum(xUp);
    upSet=cell(nXUp,1);
    
    upNxt=nxtLr(xUp)==(top-1);
    upAcro=acroLrs(xUp);
    upSet(upNxt)=upAcro(upNxt);
    
    upPrv=prvLr(xUp)==(top-1);
    upBasi=basiLrs(xUp);
    upSet(upPrv)=upBasi(upPrv);
    
    isInt=false(nXUp,1);
    for i=1:nXUp
        isInt(i)=numel(upSet{i})>0;
    end
    xUp(xUp)=isInt;
    nXUp=sum(xUp);
    
    %at bottom:
    xDn=inLayer==bot & (nxtLr==(bot+1) | prvLr==(bot+1));
    nXDn=sum(xDn);
    dnSet=cell(nXDn,1);
    
    dnNxt=nxtLr(xDn)==(bot+1);
    dnAcro=acroLrs(xDn);
    dnSet(dnNxt)=dnAcro(dnNxt);
    
    dnPrv=prvLr(xDn)==(bot+1);
    dnBasi=basiLrs(xDn);
    dnSet(dnPrv)=dnBasi(dnPrv);
    
    isSub=inLayer<=bot & inLayer>=top;
    %count minimum degrees of freedom present in each set
    
    nDofUp=redDOFs(upSet,parents,isSub,xUp,nXUp);
    nDofDn=redDOFs(dnSet,parents,isSub,xDn,nXDn);
    nDOF=nDofUp+nDofDn;
    
end