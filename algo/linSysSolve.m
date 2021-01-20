function eqs=linSysSolve(iLayer,layerEqs,prob,nLayers)

keyboard

%make all equations 'zeros' by adding depvar to vars
%then order variables in descending order of importance
%then do checks on any 'stragglers', that will ruin the solution if
%included & remove
%finally solve.

    sol=formSys(layerEqs,cat(1,layerEqs(:).kLayer),iLayer);


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
                    %in this case, should be able to use equations to cull
                    %the outer-most psiL, psiXBar
                    
                    %how to organise? this whole last part needs to be
                    %reorganised...
                    
                    
                    
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


end