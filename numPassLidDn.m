function lidEq=numPassLidDn(iSeg,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,lidEq,Kx,b2,c1,c2,c5,termed,parents,inLayer)

    par=parents(iSeg);
    
    sibs=parents==par;
    nSibs=sum(sibs); 
    
    iSibs=find(sibs);
    vars=cell(1,nSibs);
    for i=1:nSibs
        vars{i}=sprintf('G1%d',iSibs(i));
    end
    massCons=formLinkEq(par,sprintf('G0%d',par),vars,Kx(sibs)'/Kx(par));
    if nSibs==1

            if ismember(massCons.depvar,lidEq.vars)
                lidEq=subsFor(lidEq,massCons.depvar,...
                            massCons.vars,massCons.coefs);
                lidEq=sumVars(lidEq);
            end
            if strcmp(sprintf('psi0%d',par),lidEq.depvar)
                lidEq.depvar=sprintf('psi1%d',iSeg);
                lidEq=numIsolate(lidEq,'psiC');
            end
            if ismember(sprintf('G1%d',iSeg),lidEq.vars) && ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                lidEq=numPassDn(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            elseif ismember(sprintf('G1%d',iSeg),lidEq.vars)
                %this happens if passing down collar condition == psiC
                lidEq=numPassDnG(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
            elseif ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                lidEq=numPassDnPsi(lidEq,iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            end
            if strcmp(sprintf('psiC'),lidEq.depvar)
                lidEq=numIsolate(lidEq,sprintf('psi0%d',iSeg));
            end
            
    elseif nSibs==2
        sibLinkI=ismember(prob.iLinks,iSibs);
        sib=setdiff(iSibs,iSeg);
        if any(termed(sibLinkI))
            
                if ismember(massCons.depvar,lidEq.vars)
                    lidEq=subsFor(lidEq,massCons.depvar,...
                                massCons.vars,massCons.coefs);
                    lidEq=sumVars(lidEq);        
                end
                
                if strcmp(sprintf('psi0%d',par),lidEq.depvar)
                    lidEq.depvar=sprintf('psi1%d',iSeg);
                    lidEq=numIsolate(lidEq,'psiC');
                end
                if ismember(sprintf('G1%d',sib),lidEq.vars)
                    lidEq=subsFor(lidEq,closeEqs(iLinkClose==sib).depvar,...
                                closeEqs(iLinkClose==sib).vars,...
                                closeEqs(iLinkClose==sib).coefs);
                    lidEq=sumVars(lidEq);
                    lidEq.vars{strcmp(lidEq.vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                    lidEq=sumVars(lidEq);
                end
                if ismember(sprintf('G1%d',iSeg),lidEq.vars) && ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                    lidEq=numPassDn(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                elseif ismember(sprintf('G1%d',iSeg),lidEq.vars)
                    lidEq=numPassDnG(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
                elseif ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                    lidEq=numPassDnPsi(lidEq,iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                end
                if strcmp(sprintf('psiC'),lidEq.depvar)
                    lidEq=numIsolate(lidEq,sprintf('psi0%d',iSeg));
                end
        else
            iExEq=ismember(iLinkExtra,iSibs);
            nTargLegs=size(extraEqs(iExEq).targTrack,1);
            if nTargLegs==2
                keyboard
                %how to use each to go both ways?
                %create extra (lid) equations?
            elseif nTargLegs==1
                massCons=formLinkEq(iSeg,sprintf('G0%d',par),...
                    {sprintf('G1%d',iSibs(1)),sprintf('G1%d',iSibs(2))},...
                    [Kx(iSibs(1))/Kx(par) Kx(iSibs(2))/Kx(par)]);
                massCons=subsFor(massCons,extraEqs(iExEq).depvar,extraEqs(iExEq).vars,...
                    extraEqs(iExEq).coefs);
                massCons=sumVars(massCons);

                if ismember(massCons.depvar,lidEq.vars)
                    lidEq=subsFor(lidEq,massCons.depvar,...
                                massCons.vars,massCons.coefs);
                    lidEq=sumVars(lidEq);
                end
                if strcmp(sprintf('psi0%d',par),lidEq.depvar)
                    lidEq.depvar=sprintf('psi1%d',iSeg);
                    lidEq=numIsolate(lidEq,'psiC');
                end
                if ismember(extraEqs(iExEq).depvar,lidEq.vars)
                    lidEq=subsFor(lidEq,extraEqs(iExEq).depvar,...
                                extraEqs(iExEq).vars,extraEqs(iExEq).coefs);
                    lidEq=sumVars(lidEq);
                end
                if ismember(sprintf('psi1%d',sib),lidEq.vars)
                    lidEq.vars{strcmp(lidEq.vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                    lidEq=sumVars(lidEq);
                end
                if ismember(sprintf('G1%d',iSeg),lidEq.vars) && ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                    lidEq=numPassDn(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                elseif ismember(sprintf('G1%d',iSeg),lidEq.vars)
                    lidEq=numPassDnG(lidEq,iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
                elseif ismember(sprintf('psi1%d',iSeg),lidEq.vars)
                    lidEq=numPassDnPsi(lidEq,iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                end    
                if strcmp(sprintf('psiC'),lidEq.depvar)
                    lidEq=numIsolate(lidEq,sprintf('psi0%d',iSeg));
                end
            end
            
        end
    else
        keyboard
        %for >2 daughters from same junction...
        
    end
    lidEq.iLink=iSeg;
end