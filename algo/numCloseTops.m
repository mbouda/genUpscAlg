function layerEqs=numCloseTops(iSeg,layerEqs,closeEqs,iLinkClose,nLayers)

    for j=1:nLayers 
        if ismember(closeEqs(iLinkClose==iSeg).depvar,layerEqs(j).vars)
            layerEqs(j)=subsFor(layerEqs(j),closeEqs(iLinkClose==iSeg).depvar,...
                closeEqs(iLinkClose==iSeg).vars,...
                closeEqs(iLinkClose==iSeg).coefs);
            layerEqs(j)=sumVars(layerEqs(j));
        end
    end
end