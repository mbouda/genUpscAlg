function [prob,sol]=layerProbSol(iLayer,collarCond,lyrArch,params,parents,inLayer)
keyboard
    prob=numProbDef(iLayer,lyrArch.acroLrs,lyrArch.nxtLr,lyrArch.basiLrs,lyrArch.prvLr,parents,inLayer);
    
    nLayers=size(prob.kLayers,1);
    layerEqs=defLayerEqs(nLayers,prob.kLayers,inLayer,params.KrS);  %function that sets them up with KrS...
                        %should move into compLayerEqs functions for easier
                        %debugging?
    
    if any(parents(prob.iLinks)==0)
        layerEqs=compClrLyrEqs(iLayer,collarCond,layerEqs,prob,...
            params.b2,params.c1,params.c2,params.c5,params.Kx,...
            inLayer,parents,nLayers);
    else
        layerEqs=compLayerEqs(iLayer,layerEqs,prob,...
            params.b2,params.c1,params.c2,params.c5,params.Kx,...
            inLayer,parents,nLayers);
    end
    
    %sometimes psiC gets reported twice: not summed where it should be!
    mLayers=size(layerEqs,1);
    for i=1:mLayers
        layerEqs(i)=sumVars(layerEqs(i));
    end
    
    isHook=cat(1,layerEqs(:).kLayer)=='H';
    if any(isHook)
        iHook=find(isHook);
        notHook=find(~isHook)';
        for i=iHook'
            layerEqs=resolveHook(i,layerEqs,notHook);
        end
        layerEqs(isHook)=[];
    end
    mLayers=size(layerEqs,1);
    
    %can also pre-eliminate kLayer==0 cases?
    
    
    elimLayers=startsWith(cat(1,{layerEqs(:).depvar}),'psiL');
    if any(elimLayers)
        
        if sum(elimLayers)>1
            eL=find(elimLayers);
            kLayer=find(cat(1,layerEqs(:).kLayer)==iLayer);
            J=cat(2,sort(eL(eL>kLayer),'descend'),eL(eL<kLayer));
            %elminiate redundants upfront:
            keepExtra=selectExtras(layerEqs(J));
            layerEqs(J(~keepExtra))=[];
            elimLayers=startsWith(cat(1,{layerEqs(:).depvar}),'psiL');
            %should be choosing outermost layer(s) to actually keep in the
            %system, no?
        end
        mLayers=size(layerEqs,1);
        J=sort(find(elimLayers),'descend');
        
        layerEqs=truncateZeroCoeffs(layerEqs,mLayers);
        jL=1:size(layerEqs,1);
        for j=J
            %add test for duplication?? seem to have eliminated it well
            %upfront
            for k=setdiff(jL,j)
                if ismember(layerEqs(j).depvar,layerEqs(k).vars)
                    layerEqs(k)=subsFor(layerEqs(k),layerEqs(j).depvar,...
                        layerEqs(j).vars,layerEqs(j).coefs);
                    layerEqs(k)=sumVars(layerEqs(k));
                    if ismember(layerEqs(k).depvar,layerEqs(k).vars)
                        layerEqs(k)=numIsolDep(layerEqs(k));
                    end
                end
            end
            layerEqs=truncateZeroCoeffs(layerEqs,mLayers);
            jL=setdiff(jL,j);
        end
        layerEqs(J)=[];
    end
    
    sol=formSys(layerEqs,cat(1,layerEqs(:).kLayer),iLayer);
%     eqs=linSysSolve(iLayer,layerEqs,prob,nLayers); 
%     sol=eqs(prob.kLayers==iLayer);

end