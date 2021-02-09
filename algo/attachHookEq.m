function hookEq=attachHookEq(hookEq,targ,prob,Kx,b2,c1,c2,c5,parents,inLayer)

    iLink=hookEq.iLink;
    hookEq=numPassDnPsi(hookEq,iLink,inLayer(iLink),...
        c1(iLink),c2(iLink),c5(iLink));
    
    
    while iLink~=targ
        dtrs=parents==iLink;
        nDtrs=sum(dtrs);
        iDtrs=find(dtrs);
        
        if nDtrs==1
            
            par=iLink;
            iLink=iDtrs;
            
            massCons=formLinkEq(par,sprintf('G0%d',par),...
                sprintf('G1%d',iLink),Kx(iLink)/Kx(par));
            hookEq=subsFor(hookEq,massCons.depvar,...
                                massCons.vars,massCons.coefs);
            hookEq.vars{strcmp(hookEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iLink);
            if iLink~=targ
                hookEq=numPassDn(hookEq,iLink,inLayer(iLink),b2(iLink),c1(iLink),c2(iLink),c5(iLink));
            end
            hookEq.iLink=iLink;
        elseif nDtrs==2
            lands=false(2,1);
            for j=1:2
                lands(j)=srchLegTarg(iDtrs(j),parents,targ,prob.terms);
            end
            
            %currently assumes that sib terms in layers above domain...
            %if fails, get psiL in here from the domain...
                %but if fails worse, get psiL from BELOW!
                %which adds unknowns/undesired psiL. 
            %so, eventually, should here stop at upper border of domain
                %add checks to the functions used above...
                %and then keep those unknowns? they should be tops in the
                %problem or other hangers... so they should end up with an
                %equation for elimination.
            
            par=iLink;
            
            sib=iDtrs(~lands);
            iLink=iDtrs(lands);
            
            if any(parents==sib)
                j=distalTipSrch(sib,parents);
            else
                j=sib;
            end
            nTerms=size(j,1);
            [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
            [closeEqs,iLinkClose,targTrack]=numCloseToSeg(par,closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
            
            %now, will form mass Cons, use closeEq to subs
            %subs into hookEq and pass down
            massCons=formLinkEq(par,sprintf('G0%d',par),...
                {sprintf('G1%d',iLink), sprintf('G1%d',sib)},...
                [Kx(iLink)/Kx(par) Kx(sib)/Kx(par)]);
            
            hookEq=subsFor(hookEq,massCons.depvar,...
                                massCons.vars,massCons.coefs);
            hookEq=subsFor(hookEq,closeEqs.depvar,...
                                closeEqs.vars,closeEqs.coefs);
            hookEq.vars{strcmp(hookEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iLink);
            hookEq.vars{strcmp(hookEq.vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iLink);
            hookEq=sumVars(hookEq);        
            if iLink~=targ
                hookEq=numPassDn(hookEq,iLink,inLayer(iLink),b2(iLink),c1(iLink),c2(iLink),c5(iLink));
            end
            hookEq.iLink=iLink;
        else
            keyboard
        end
        
    end

end