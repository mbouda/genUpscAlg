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
        iDtrs=find(dtrs);
        sib=setdiff(iDtrs,hookEq.iLink);
        
        if hookEq.iLink~=top
            j=distalTipSrch(sib,parents);
            nTerms=size(j,1);
            [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
            [closeEqs,iLinkClose,targTrack]=numCloseToSeg(parents(sib),closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
            massCons=formLinkEq(par,sprintf('G1%d',hookEq.iLink),...
                    {sprintf('G0%d',par) sprintf('G1%d',sib)},...
                    [Kx(par)/Kx(hookEq.iLink) Kx(sib)/Kx(hookEq.iLink)]);
            massCons=subsFor(massCons,closeEqs.depvar,closeEqs.vars,closeEqs.coefs);
            hookEq=subsFor(hookEq,massCons.depvar,...
                massCons.vars,massCons.coefs);
            hookEq.vars{ismember(hookEq.vars,sprintf('psi1%d',sib))}=sprintf(sprintf('psi1%d',hookEq.iLink));
            hookEq=sumVars(hookEq);
            hookEq=numIsolDep(hookEq);
            hookEq=numPassUpG(hookEq,par,inLayer(par),b2(par),c1(par),c2(par));
            hookEq.iLink=par;
        end
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
            iDtrs=find(dtrs);
            sib=setdiff(iDtrs,hookEq.iLink);
            j=distalTipSrch(sib,parents);
            nTerms=size(j,1);
            [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
            [closeEqs,iLinkClose,targTrack]=numCloseToSeg(parents(sib),closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
            
            massCons=formLinkEq(par,sprintf('G1%d',hookEq.iLink),...
                {sprintf('G0%d',par) sprintf('G1%d',sib)},...
                [Kx(par)/Kx(hookEq.iLink) Kx(sib)/Kx(hookEq.iLink)]);
            massCons=subsFor(massCons,closeEqs.depvar,closeEqs.vars,closeEqs.coefs);
            hookEq=subsFor(hookEq,massCons.depvar,...
                massCons.vars,massCons.coefs);
            hookEq.vars{ismember(hookEq.vars,sprintf('psi1%d',hookEq.iLink))}=sprintf(sprintf('psi0%d',par));
            hookEq.vars{ismember(hookEq.vars,sprintf('psi1%d',sib))}=sprintf(sprintf('psi0%d',par));
            hookEq=sumVars(hookEq);
            
            hookEq=numPassUp(hookEq,par,inLayer(par),b2(par),c1(par),c2(par),c5(par));
            hookEq.iLink=par;
            
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