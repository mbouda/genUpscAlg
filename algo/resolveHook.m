function layerEqs=resolveHook(iHook,layerEqs,notHook)

    for i=notHook
        if ismember(layerEqs(iHook).depvar,layerEqs(i).vars)
            layerEqs(i)=subsFor(layerEqs(i),layerEqs(iHook).depvar,...
                layerEqs(iHook).vars,layerEqs(iHook).coefs);
            layerEqs(i)=sumVars(layerEqs(i));
        end
    end
end