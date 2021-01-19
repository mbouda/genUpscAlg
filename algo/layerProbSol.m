function [prob,sol]=layerProbSol(iLayer,collarCond,lyrArch,params,parents,inLayer)

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
        end
        J=sort(find(elimLayers),'descend');
        
        jL=1:size(layerEqs,1);
        for j=J
            %add test for duplication??
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
            if any(layerEqs(j).coefs==0) || any(abs(layerEqs(j).coefs)==1) || ...
                    any(abs(layerEqs(j).coefs)<(1e-16*max(abs(layerEqs(j).coefs)))) 
            else
                warning('duplicate equation? zero/one coeffs','degenCoefs')
                %shoulad actually likely just remove these coefficients
                %before this loop...
                %use: layerEqs=truncateZeroCoeffs(layerEqs,nLayers);
            end
            jL=setdiff(jL,j);
        end
        layerEqs(J)=[];
    end
    
    
    
    
    [elimVars,nVars]=countSysVars(layerEqs(:),prob);
    mLayers=size(layerEqs,1);
    
    if nVars==(mLayers-1)
        eqs=solveSysFor(iLayer,prob.kLayers,prob.kLayers,layerEqs,elimVars,nVars);
    elseif nVars<(nLayers-1)

        layAbs=cell(nLayers,1);
        for i=1:nLayers
            psiLI=discoverIndices(layerEqs(i).vars,'psiL');
            layAbs{i}=setdiff(prob.kLayers,psiLI);
        end
        
        %test for odd ones out, can eliminate that equation...
        layerSet=prob.kLayers;
        nEq=nLayers-1;
        j=1;
        while nEq>nVars && j<=nEq
            
            isAbs=false(nLayers,1);
            for k=1:nLayers
                isAbs(k)=ismember(layerSet(j),layAbs{k});
            end
            
            if sum(isAbs)>=(nVars+1)
                layerSet(j)=[];
                nEq=nEq-1;
            else
                j=j+1;
            end
        end
        if ~ismember(iLayer,layerSet)
            layerSet=cat(1,layerSet(layerSet<iLayer),iLayer,layerSet(layerSet>iLayer));
            nEq=nEq+1;
        end
        if nVars==nEq
            eqs=solveSysFor(iLayer,layerSet,prob.kLayers,layerEqs,elimVars,nVars);
        else
            %could still try to minimise variables in resulting equation
            %here..
            isInSet=false(nLayers,1);
            isInSet(prob.kLayers==iLayer)=true;
            down=find(isInSet)+1;
            up=find(isInSet)-1;
            while sum(isInSet)<(nVars+1) && sum(isInSet)<nLayers
                if down<=nLayers
                    isInSet(down)=true;
                    down=down+1;
                end
                
                if sum(isInSet)<(nVars+1) && up>0
                    isInSet(up)=true;
                    up=up-1;
                end
            end
            layerSet=prob.kLayers(isInSet);
            eqs=solveSysFor(iLayer,layerSet,prob.kLayers,layerEqs,elimVars,nVars);
        end
    else
        
        layerEqs=truncateZeroCoeffs(layerEqs,nLayers);
        [elimVars,nVars]=countSysVars(layerEqs(:),prob);
        if nVars>(nLayers-1)
            %underdetermined
            %for lupinus, 42-day, can solve if drop layer(s)
            %insert code testing subsets of layers
            jLayer=find(prob.kLayers==iLayer);
            jTop=jLayer;
            jBot=jLayer;

            [elimVars,nVars]=countSysVars(layerEqs(jTop:jBot),prob);
            while nVars>jBot-jTop && (jTop > 1 || jBot<nLayers)
                %probably house this in separate function:
                kBot=min(jBot+1,nLayers);
                kTop=max(jTop-1,1);

                [elimVars,nVars]=countSysVars(layerEqs(jTop:kBot),prob);
                if nVars<=kBot-jTop
                    jBot=kBot;
                else
                    [elimVars,nVars]=countSysVars(layerEqs(kTop:jBot),prob);
                    if nVars<=jBot-kTop
                        jTop=kTop;
                    else
                        [elimVars,nVars]=countSysVars(layerEqs(kTop:kBot),prob);
                        jBot=kBot;
                        jTop=kTop;
                    end
                end
                %this algorithm steps out symmetrically from the iLayer
                %fails to account for cases where system is solvable by moving
                %asymmetrically up or down... would have to go through further
                %algorithm to just scan all (contiguous?) possibilities that
                %overlap the iLayer
            end
            if nVars<=jBot-jTop
                eqs=solveSysFor(iLayer,(jTop:jBot)',prob.kLayers,layerEqs(jTop:jBot),elimVars,nVars);
            else
                if nVars<=size(layerEqs,1)
                    keyboard
                else
                    keyboard
                    %still underdetermined: algorithm failed
                    %usually down to a segment ordering issue...
                end
            end
        elseif nVars==(nLayers-1)
            eqs=solveSysFor(iLayer,prob.kLayers,prob.kLayers,layerEqs(:),elimVars,nVars);
        else
            keyboard
            %now overdetermined....
        end
    end
    
    sol=eqs(prob.kLayers==iLayer);

end