function [hookEqs,nHL]=connectHangingTarg(hangLinks,closeEqs,iLinkClose,prob,parents,b2,c1,c2,c5,Kx,inLayer)

        nHL=size(hangLinks,1);
        
        isConn=false(nHL,1);
        tops=zeros(nHL,1);
        hits=zeros(nHL,1);
        pars=parents(hangLinks);
        
        %alg based on bad ass'n: max(pars) may already be connected, no?
        while any(~isConn)
            
            [~,imp]=max(pars);
            
            [touches,hit]=srchTarg(pars(imp),parents,cat(1,setdiff(prob.targ,(hangLinks(imp))),hangLinks(~imp)),cat(1,prob.bots,prob.terms,hangLinks(imp)));
            if touches
                isConn(imp)=true;
                hits(imp)=hit;
                tops(imp)=pars(imp);  
                pars(imp)=0;
            else
                dtrs(imp)=pars(imp);
                pars(imp)=parents(pars(imp));
            end
        end
        
        %tops end up being the pars, who's daughters connect down to targets.
        
        hookEqs=struct('depvar',cell(nHL,1),'vars',cell(nHL,1),'coefs',cell(nHL,1),'iLink',cell(nHL,1));
        
        for i=1:nHL
            par=tops(i);
            sib=setdiff(find(parents==par),dtrs(i));
            
            iClEq=(iLinkClose==hangLinks(i));
            hookEqs(i)=formHookEq(closeEqs(iClEq),hangLinks(i),dtrs(i),prob,Kx,b2,c1,c2,c5,parents,inLayer);
            hookEqs(i).vars{ismember(hookEqs(i).vars,sprintf('psi1%d',dtrs(i)))}=sprintf(sprintf('psi1%d',sib));
            hookEqs(i).iLink=sib;
            
            hookEqs(i)=attachHookEq(hookEqs(i),hits(i),prob,Kx,b2,c1,c2,c5,parents,inLayer);
        end
        
        
end