function [layerEqs,prob,nLayers]=numConnUD(iSeg,layerEqs,closeEqs,iLinkClose,extraEqs,iLinkExtra,nLayers,prob,termed,Kx,parents)

    par=parents(iSeg);

    sibs=parents==par;
    nSibs=sum(sibs);  %including self
    iSibs=find(sibs);

    if nSibs==1
        for j=1:nLayers
            if ismember(sprintf('G0%d',par),layerEqs(j).vars)
                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                layerEqs(j)=subsFor(layerEqs(j),sprintf('G0%d',par),...
                            {sprintf('G1%d',iSeg)},...
                            Kx(iSeg)/Kx(par));
                layerEqs(j)=sumVars(layerEqs(j));
            end
        end

    elseif nSibs==2

        sib=setdiff(iSibs,iSeg);
        massCons=formLinkEq(iSeg,sprintf('G0%d',par),...
                {sprintf('G1%d',iSibs(1)),sprintf('G1%d',iSibs(2))},...
                [Kx(iSibs(1))/Kx(par) Kx(iSibs(2))/Kx(par)]);   
        if any(termed(ismember(prob.iLinks,iSibs))) %should never be case that all are closed.
            isSib=ismember(iLinkClose,sib);
                %actually, should be just in one position and sib is scalar
                %so maybe can do iLinkClose==sib here for efficiency?
            for j=1:nLayers
                isJ=prob.kLayers==prob.kLayers(j);
                if ismember(massCons.depvar,layerEqs(isJ).vars)
                    layerEqs(isJ)=subsFor(layerEqs(isJ),massCons.depvar,...
                                massCons.vars,massCons.coefs);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
                
                if ismember(closeEqs(isSib).depvar,layerEqs(isJ).vars)
                    layerEqs(isJ)=subsFor(layerEqs(isJ),closeEqs(isSib).depvar,...
                                closeEqs(isSib).vars,...
                                closeEqs(isSib).coefs);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end 
                if ismember(sprintf('psi1%d',sib),layerEqs(isJ).vars) 
                    layerEqs(isJ).vars{strcmp(layerEqs(isJ).vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
                if ismember(sprintf('psi0%d',par),layerEqs(isJ).vars)
                    layerEqs(isJ).vars{strcmp(layerEqs(isJ).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
            end
        else
            if ismember(sib,prob.targ)
                keyboard
                %need extra code for if both are targets...
                %may just need to re-solve extraEq for iSeg once out of the two sibs
            end
            isSib=ismember(iLinkExtra,sib);
            for j=1:nLayers
                isJ=prob.kLayers==prob.kLayers(j);
                if ismember(massCons.depvar,layerEqs(isJ).vars)
                    layerEqs(isJ)=subsFor(layerEqs(isJ),massCons.depvar,...
                                massCons.vars,massCons.coefs);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end 
                if ismember(sprintf('psi0%d',par),layerEqs(isJ).vars)
                    layerEqs(isJ).vars{strcmp(layerEqs(isJ).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
                if ismember(extraEqs(isSib).depvar,layerEqs(isJ).vars)
                    layerEqs(isJ)=subsFor(layerEqs(isJ),extraEqs(isSib).depvar,...
                                extraEqs(isSib).vars,...
                                extraEqs(isSib).coefs);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
                if ismember(sprintf('psi1%d',sib),layerEqs(isJ).vars)  %does this ever occur?
                    layerEqs(isJ).vars{strcmp(layerEqs(isJ).vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                    layerEqs(isJ)=sumVars(layerEqs(isJ));
                end
            end
            [layerEqs,prob,nLayers]=addLayerEq(extraEqs(isSib),layerEqs,...
                    prob,nLayers,massCons,iSeg,par,'ud',...
                    [],[],[],[],[]);
        end
    else
        keyboard
        %for multiple daughters from same junction...

    end

end