function [layerEqs,prob,nLayers]=addLayerEq(extraEq,layerEqs,prob,nLayers,massCons,iSeg,par,dir,inLayer,b2,c1,c2,c5)

    extraLayers=discoverIndices(extraEq.vars,'psiL');
    extraLayers=setdiff(extraLayers,prob.kLayers);
    if any(extraLayers)
        newLayer=extraLayers(1);
        newEq=numIsolate(extraEq,sprintf('psiL%d',newLayer));
        newEq=rmfield(newEq,'helperEqs');
        newEq=rmfield(newEq,'targTrack');

        switch dir
            case 'up'
                massCons=numIsolate(massCons,sprintf('G1%d',extraEq.iLink));
                newEq=subsFor(newEq,massCons.depvar,massCons.vars,massCons.coefs);
                newEq=sumVars(newEq);
                if ismember(newEq.depvar,newEq.vars)
                    newEq=numIsolDep(newEq);
                end
                isSib=cat(1,extraEq.helperEqs(:).iLink)~=extraEq.iLink;
                newEq=subsFor(newEq,extraEq.helperEqs(isSib).depvar,...
                    extraEq.helperEqs(isSib).vars,extraEq.helperEqs(isSib).coefs);  
                newEq=sumVars(newEq);
                if ismember(newEq.depvar,newEq.vars)
                    newEq=numIsolDep(newEq);
                end
                newEq=numPassUp(newEq,iSeg,inLayer,b2,c1,c2,c5);
            case 'dn'
                isSib=cat(1,extraEq.helperEqs(:).iLink)~=iSeg;
                newEq=subsFor(newEq,extraEq.helperEqs(isSib).depvar,...
                    extraEq.helperEqs(isSib).vars,extraEq.helperEqs(isSib).coefs);  
                newEq=sumVars(newEq);
                if ismember(newEq.depvar,newEq.vars)
                    if newEq.coefs(ismember(newEq.vars,newEq.depvar))==1
                        
                        %in this case, do not use the extraEq here, as it
                        %will return NaNs upun numIsolDep
                        
                        %or, try: eliminating the depvar from vars, make depvar 0
                        newEq.vars(ismember(newEq.vars,newEq.depvar))=[];
                        newEq.coefs(ismember(newEq.vars,newEq.depvar))=[];
                        
                        newEq.depvar='0';
                        newLayer=0;
                        
                        %or select a different layer off the list
                        
                        
                        
                        %return
                        
                        
                    else
                        newEq=numIsolDep(newEq);
                    end
                end
                newEq.vars{ismember(newEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                newEq=sumVars(newEq);
                newEq=numPassDn(newEq,iSeg,inLayer,b2,c1,c2,c5);
            case 'ud'
                isSib=cat(1,extraEq.helperEqs(:).iLink)~=iSeg;
                newEq=subsFor(newEq,extraEq.helperEqs(isSib).depvar,...
                    extraEq.helperEqs(isSib).vars,extraEq.helperEqs(isSib).coefs);  
                newEq=sumVars(newEq);
                newEq=numIsolDep(newEq);

                newEq.vars{ismember(newEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                newEq=sumVars(newEq);
        end
        if ismember(newEq.depvar,newEq.vars)
            newEq=numIsolDep(newEq);
        end
            %ADD newEq to layerEqs
        newEq=rmfield(newEq,'iLink');
        newEq.kLayer=newLayer;

        layerEqs=cat(1,layerEqs,newEq);

        prob.kLayers=cat(1,prob.kLayers,newLayer);
        nLayers=nLayers+1;
    elseif numel(extraEq.targTrack)>1
        newEq=extraEq;
        newEq=rmfield(newEq,'helperEqs');
        newEq=rmfield(newEq,'targTrack');
        newEq=rmfield(newEq,'iLink');

        isSib=cat(1,extraEq.helperEqs(:).iLink)==iSeg;
        newEq=subsFor(newEq,extraEq.helperEqs(isSib).depvar,...
            extraEq.helperEqs(isSib).vars,extraEq.helperEqs(isSib).coefs);
        newEq=sumVars(newEq);
        newEq.vars{ismember(newEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
        newEq=sumVars(newEq);
        newEq=numPassDnPsi(newEq,iSeg,inLayer,c1,c2,c5);

        newEq.vars=cat(2,newEq.vars,newEq.depvar);
        newEq.coefs=cat(2,newEq.coefs,-1);
        newEq.depvar='0';

        newEq=subsFor(newEq,extraEq.helperEqs(~isSib).depvar,...
            extraEq.helperEqs(~isSib).vars,extraEq.helperEqs(~isSib).coefs);
        newEq.vars{ismember(newEq.vars,sprintf('psi0%d',par))}=sprintf('psi1%d',...
            extraEq.helperEqs(~isSib).iLink);
        newEq=sumVars(newEq);

        newEq=truncateZeroCoeffs(newEq,1);

        newEq.kLayer=0;
        layerEqs=cat(1,layerEqs,newEq);
        prob.kLayers=cat(1,prob.kLayers,0);
        nLayers=nLayers+1;
    end
end