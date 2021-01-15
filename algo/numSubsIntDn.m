function [layerEqs,prob,nLayers]=numSubsIntDn(iSeg,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,layerEqs,nLayers,Kx,b2,c1,c2,c5,termed,parents,inLayer)

    iLayer=prob.kLayers==inLayer(iSeg);
    if ismember(sprintf('psiX%d',iSeg),layerEqs(iLayer).vars)
        layerEqs(iLayer)=subsIntoAvg(layerEqs(iLayer),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg));
    end
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
        for j=1:nLayers
            if ismember(massCons.depvar,layerEqs(j).vars)
                layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,...
                            massCons.vars,massCons.coefs);
                layerEqs(j)=sumVars(layerEqs(j));
            end
            if ismember(sprintf('psi0%d',par),layerEqs(j).vars)
                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                layerEqs(j)=sumVars(layerEqs(j));
            end
            if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                layerEqs(j)=numPassDn(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            elseif ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                %this happens if passing down collar condition == psiC
                layerEqs(j)=numPassDnG(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
            elseif ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                layerEqs(j)=numPassDnPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            end
        end
    elseif nSibs==2
        sibLinkI=ismember(prob.iLinks,iSibs);
        sib=setdiff(iSibs,iSeg);
        if any(termed(sibLinkI))
            for j=1:nLayers
                if ismember(massCons.depvar,layerEqs(j).vars)
                    layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,...
                                massCons.vars,massCons.coefs);
                    layerEqs(j)=sumVars(layerEqs(j));        
                end
                if ismember(sprintf('psi0%d',par),layerEqs(j).vars)
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                end
                if ismember(sprintf('G1%d',sib),layerEqs(j).vars)
                    layerEqs(j)=subsFor(layerEqs(j),closeEqs(iLinkClose==sib).depvar,...
                                closeEqs(iLinkClose==sib).vars,...
                                closeEqs(iLinkClose==sib).coefs);
                    layerEqs(j)=sumVars(layerEqs(j));
                end
                if ismember(sprintf('psi1%d',sib),layerEqs(j).vars)
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                    layerEqs(j)=sumVars(layerEqs(j));
                end
                if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                    layerEqs(j)=numPassDn(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                elseif ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                    layerEqs(j)=numPassDnPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                elseif ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                    layerEqs(j)=numPassDnG(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
                end
            end
        else
            iExEq=ismember(iLinkExtra,iSibs);
            nTargLegs=size(extraEqs(iExEq).targTrack,1);
            if nTargLegs==2
                keyboard
                %how to use each to go both ways?
            elseif nTargLegs==1
                massCons=formLinkEq(iSeg,sprintf('G0%d',par),...
                    {sprintf('G1%d',iSibs(1)),sprintf('G1%d',iSibs(2))},...
                    [Kx(iSibs(1))/Kx(par) Kx(iSibs(2))/Kx(par)]);
                massCons=subsFor(massCons,extraEqs(iExEq).depvar,extraEqs(iExEq).vars,...
                    extraEqs(iExEq).coefs);
                massCons=sumVars(massCons);

                for j=1:nLayers
                    if ismember(massCons.depvar,layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,...
                                    massCons.vars,massCons.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('psi0%d',par),layerEqs(j).vars) %apparently fails when passing down psiC
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(extraEqs(iExEq).depvar,layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),extraEqs(iExEq).depvar,...
                                    extraEqs(iExEq).vars,extraEqs(iExEq).coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('psi1%d',sib),layerEqs(j).vars)
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassDn(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    elseif ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassDnPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    elseif ismember(sprintf('G1%d',iSeg),layerEqs(j).vars) 
                        layerEqs(j)=numPassDnG(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg));
                    end
                    if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                        layerEqs(j)=numIsolDep(layerEqs(j));
                    end
                end
                [layerEqs,prob,nLayers]=addLayerEq(extraEqs(iExEq),layerEqs,...
                    prob,nLayers,massCons,iSeg,par,'dn',...
                    inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            end
            
        end
    else
        keyboard
        %for >2 daughters from same junction...
        
    end

end