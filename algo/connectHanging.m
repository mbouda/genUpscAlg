function hookEqs=connectHanging(hangLinks,closeEqs,iLinkClose,prob,parents,b2,c1,c2,c5,Kx,inLayer)
keyboard
        nHL=size(hangLinks,1);
        
        isConn=false(nHL,1);
        tops=zeros(nHL,1);
        hits=zeros(nHL,1);
        pars=parents(hangLinks);
        
        %alg based on ass'n: max(pars) may already be connected, no?
        while any(~isConn)
            
            [~,imp]=max(pars);
            
            [touches,hit]=srchTarg(pars(imp),parents,cat(1,prob.targ,hangLinks(~imp)),cat(1,prob.bots,prob.terms,hangLinks(imp)));
            if touches
                isConn(imp)=true;
                hits(imp)=hit;
                tops(imp)=pars(imp);  
                pars(imp)=0;
            else
                pars(imp)=parents(pars(imp));
            end
        end
        
        %tops end up being the links whos SIBs connect down to targets.
        
        
        hookEqs=struct('depvar',cell(nHL,1),'vars',cell(nHL,1),'coefs',cell(nHL,1),'iLink',cell(nHL,1));
        
        for i=1:nHL
            par=parents(tops(i));
            sib=setdiff(find(parents==par),tops(i));
            
            j=distalTipSrch(tops(i),parents);
            k=distalTipSrch(hangLinks(i),parents);
            j=setdiff(j,k);
            if any(j)
                nTerms=size(j,1);
                [locCloseEqs,jLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
            end
            
            
            iClEq=(iLinkClose==hangLinks(i));
            hookEqs(i)=formHookEq(closeEqs(iClEq),hangLinks(i),tops(i),prob,Kx,b2,c1,c2,c5,parents,inLayer);
            
            hookEqs(i).vars{ismember(hookEqs(i).vars,sprintf('psi1%d',tops(i)))}=sprintf(sprintf('psi1%d',sib));
            hookEqs(i).iLink=sib;
            hookEqs(i)=attachHookEq(hookEqs(i),hits(i),prob,Kx,b2,c1,c2,c5,parents,inLayer);
        end
        
        hookEqs=rmfield(hookEqs,'iLink');
        [hookEqs(:).kLayer]=deal('H');
end