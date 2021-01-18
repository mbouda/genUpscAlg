function [extraEq,iLinkEx]=compExtraEq(iLink,prob,b2,c1,c2,c5,Kx,inLayer,parents)

    j=distalTipSrch(iLink,parents);
    nTerms=size(j,1);
    [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
    [closeEqs,iLinkClose,targTrack]=numCloseToSeg(iLink,closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
    
    if size(iLinkClose,1)>2
        keyboard
    end
    
    for i=1:2 
        closeEqs(i).vars{strcmp(closeEqs(i).vars,sprintf('psi1%d',iLinkClose(i)))}=sprintf('psi0%d',iLink);
    end
    
    hasTargs=cellfun(@(x)any(x),{targTrack(:).list})';
    
    if ~any(hasTargs)
        eqSub=numIsolate(closeEqs(1),sprintf('psi0%d',iLink));
        extraEq=subsFor(closeEqs(2),eqSub.depvar,eqSub.vars,eqSub.coefs);
        extraEq=sumVars(extraEq);
        extraEq.helperEqs=closeEqs;
        iLinkEx=iLinkClose(2);
        extraEq.targTrack=[];
    elseif all(hasTargs)
        eqSub=numIsolate(closeEqs(1),sprintf('psi0%d',iLink));
        extraEq=subsFor(closeEqs(2),eqSub.depvar,eqSub.vars,eqSub.coefs);
        extraEq=sumVars(extraEq);
        extraEq.helperEqs=closeEqs;
        iLinkEx=iLinkClose(2);
        extraEq.targTrack=targTrack;
    else
        targHead=ismember(iLinkClose,targTrack(hasTargs).head);
        eqSub=numIsolate(closeEqs(targHead),sprintf('psi0%d',iLink));
        extraEq=subsFor(closeEqs(~targHead),eqSub.depvar,eqSub.vars,eqSub.coefs);
        extraEq=sumVars(extraEq);
        extraEq.helperEqs=closeEqs;
        iLinkEx=iLinkClose(~targHead);
        extraEq.targTrack=targTrack(targHead);
    end
end