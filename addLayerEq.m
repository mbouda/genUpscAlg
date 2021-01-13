function [layerEqs,prob,nLayers]=addLayerEq(extraEq,layerEqs,prob,nLayers,massCons,iSeg,par,dir,inLayer,b2,c1,c2,c5)
        
%for multiple added eqs in series, appears to result in
%degenerate coefficients
%pretty much regardless of how it is done...


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
                        newEq=numIsolDep(newEq);
                        
                        newEq=numPassUp(newEq,iSeg,inLayer,b2,c1,c2,c5);
                    case 'dn'
                        isSib=cat(1,extraEq.helperEqs(:).iLink)~=iSeg;
                        newEq=subsFor(newEq,extraEq.helperEqs(isSib).depvar,...
                            extraEq.helperEqs(isSib).vars,extraEq.helperEqs(isSib).coefs);  
                        newEq=sumVars(newEq);
                        newEq=numIsolDep(newEq);
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
            end
end