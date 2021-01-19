function hookEq=formHookEq(eqIn,hanger,top,prob,Kx,b2,c1,c2,c5,parents,inLayer)

    hookEq=numIsolate(eqIn,sprintf('psi1%d',hanger));
    
    par=parents(hookEq.iLink);
    dtrs=parents==par;
    nDtrs=sum(dtrs);
       
    if nDtrs==1
        hookEq=subsFor(hookEq,sprintf('G1%d',hanger),...
            sprintf('G0%d',par),...
            Kx(par)/Kx(hanger));
        hookEq=numPassUpG(hookEq,par,inLayer(par),b2(par),c1(par),c2(par));
        hookEq.iLink=par;
    elseif nDtrs==2
        keyboard
    else
        keyboard
    end
    
    
    while hookEq.iLink~=top
    
        par=parents(hookEq.iLink);
        dtrs=parents==par;
        nDtrs=sum(dtrs);
        
        if nDtrs==1
            hookEq=subsFor(hookEq,sprintf('G1%d',hookEq.iLink),...
                sprintf('G0%d',par),...
                Kx(par)/Kx(hookEq.iLink));
            hookEq.vars{ismember(hookEq.vars,sprintf('psi1%d',hookEq.iLink))}=sprintf(sprintf('psi0%d',par));
            hookEq=numPassUp(hookEq,par,inLayer(par),b2(par),c1(par),c2(par),c5(par));
            hookEq.iLink=par;
        elseif nDtrs==2
            keyboard
            
            iDtrs=find(dtrs);
            sib=setdiff(iDtrs,hookEq.iLink);
            
            j=distalTipSrch(sib,parents);
            
            nTerms=size(j,1);
            [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
            [closeEqs,iLinkClose,targTrack]=numCloseToSeg(sib,closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
            %should result in close equation up to the sib (may throw error at end... 
                %looks like it's designed for parent, but don't want to duplicate calculations...)
            
            
        else
            keyboard
            %deal with supernumerary juncitons here
        end
        
    end
    
    j=distalTipSrch(top,parents);
            
    nTerms=size(j,1);
    [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
    [closeEqs,iLinkClose,targTrack]=numCloseToSeg(parents(top),closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
    
    hookEq=subsFor(hookEq,closeEqs.depvar,closeEqs.vars,closeEqs.coefs);
    hookEq=sumVars(hookEq);
    
end