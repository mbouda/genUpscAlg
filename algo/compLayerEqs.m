function layerEqs=compLayerEqs(iLayer,layerEqs,prob,b2,c1,c2,c5,Kx,inLayer,parents,nLayers)
keyboard
    topTerms=intersect(prob.terms,prob.tops);
    topTermSibs=[];
    if any(topTerms)
        prob.tops(prob.tops==topTerms)=[];
        for j=topTerms'
            topTermSibs=cat(2,topTermSibs,setdiff(find(parents==parents(j)),j));
        end
    end

    j=prob.tops;
    topLK=ismember(prob.kLayers,inLayer(j));
    layerEqs(topLK)=subsTops(j,iLayer,layerEqs(topLK),c1(j),c2(j),prob.kLayers(topLK),inLayer(j),numel(j));
    for j=find(topLK)'
        layerEqs(j)=sumVars(layerEqs(j));
    end

    j=prob.bots;
    botLK=ismember(prob.kLayers,inLayer(j));
    layerEqs(botLK)=subsBots(j,layerEqs(botLK),c1(j),c2(j),prob.kLayers(botLK),inLayer(j),numel(j));
    for j=find(botLK)'
        layerEqs(j)=sumVars(layerEqs(j));
    end


    j=prob.terms;
    nTerms=size(j,1);
    [closeEqs,iLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
    termed=false(size(prob.iLinks));
    termed(ismember(prob.iLinks,j))=true;

    termLK=ismember(prob.kLayers,inLayer(j));
    layerEqs(termLK)=subsBots(j,layerEqs(termLK),c1(j),c2(j),prob.kLayers(termLK),inLayer(j),numel(j));
    for j=find(termLK)'
        layerEqs(j)=sumVars(layerEqs(j));
    end

    [closeEqs,layerEqs,termed,iLinkClose]=numCloseInts(prob,closeEqs,iLinkClose,layerEqs,nLayers,Kx,b2,c1,c2,c5,...
                                     parents,inLayer,termed,iLayer);
    
    oJuncPars=getOpenJuncs(prob.iLinks,parents,termed);
    if any(oJuncPars)
        nOJ=size(oJuncPars,1);
        
%         if nOJ>1
%             %for efficiency, can relate junctions to their junction
%             %parent, so that can save on repeating calculations up the
%             %network by taking an acropetal extraEq and propagating
%             %it basipetally...
%             
%             keyboard
%             %From each unclosed junction, go up parents until meet another or a top; 
%                 %record relations and which leg arrive at.
%         end

        extraEqs=struct('iLink',cell(nOJ,1),'depvar',cell(nOJ,1),'vars',cell(nOJ,1),...
                    'coefs',cell(nOJ,1),'helperEqs',cell(nOJ,1),'targTrack',cell(nOJ,1));
        iLinkExtra=zeros(nOJ,1);
        for j=1:nOJ
            %in the compExtraEq function, need to track which (if any) sibs
            %have targs below them, and which ones(?)
            [extraEqs(j),iLinkExtra(j)]=compExtraEq(oJuncPars(j),prob,b2,c1,c2,c5,Kx,inLayer,parents);
        end
    else
        extraEqs=[];
        iLinkExtra=[];
        nOJ=0;
    end
           
    %identify hanging elements
    allTops=union(prob.tops,topTerms);
    hangLinks=identifyHanging(allTops,prob,parents);
    if any(hangLinks)
            hangClosed=ismember(hangLinks,iLinkClose);
            hookEqs=connectHanging(hangLinks(hangClosed),closeEqs,iLinkClose,prob,parents,b2,c1,c2,c5,Kx,inLayer);
            %construct minimal network that connects each hanging Link to 1(+)
            %targ or one other hanging link...
            %need one equation for each hanger
            
        if any(~hangClosed)
            hookUnclEqs=connectHangingUnclosed(hangLinks(~hangClosed),prob,parents,b2,c1,c2,c5,Kx,inLayer);
            hookEqs=cat(1,hookEqs,hookUnclEqs);
        else
            
        end
    else
        hookEqs=[];
    end
    
    % this assumes in a strong way that acropetal is down, basipetal is up
        %fails for pisum sativum... need to traverse actual network.. 

    %all others are upLinks?
    downLinks=idDownLinks(prob.targ,allTops',parents);
    downInts=intersect(downLinks,prob.ints)';
    upInts=setdiff(setdiff(prob.ints,downInts),prob.iLinks(termed))';
    
    
    %downInts=intersect(prob.ints(inLayer(prob.ints)<iLayer),prob.iLinks(~termed))';
    %upInts=flipud(intersect(prob.ints(inLayer(prob.ints)>=iLayer),prob.iLinks(~termed)))';
    
    for j=1:nOJ
        nLegsTarg=numel(extraEqs(j).targTrack);
        if nLegsTarg==1
            downExtra=ismember(downInts,iLinkExtra(j));
            downExtra=runDownExtras(iLinkExtra(j),downExtra,downInts,parents);
            topExtra=ismember(prob.tops,iLinkExtra(j));
            downExtra=runDownExtras(iLinkExtra(j),downExtra,downInts,parents);
            if any(topExtra)
                upInts=cat(2,upInts,prob.tops(topExtra));
            end
            upInts=cat(2,upInts,downInts(downExtra));
            downInts=downInts(~downExtra);
            sib=setdiff(cat(1,extraEqs(j).helperEqs(:).iLink),extraEqs(j).iLink);
            if ~ismember(sib,downInts)
                downInts=cat(2,downInts,sib);
            end
%         elseif nLegsTarg==2
%             keyboard
            
            %will need to be both up AND down?
            %take both siblings
            %do not take off down list
            %add to up list? (along with relevant descendants)
            %will need to subs on way down! but maybe doing upward is
            %duplication of extraEq work... since we know it goes to targ
            %anyway...
            
            %seems like only need to add to up list descendants that do NOT
            %terminate in targ (?), i.e. just do the first part of the
            %conditional
            
        end %if ==0, no need to execute anything
    end
    downInts=cat(2,downInts,topTermSibs);
    
    upInts=sort(upInts,'descend');
    downInts=sort(downInts,'ascend');

    if size(upInts,1)~=1 %for some reason, needed now...
        upInts=upInts';
    end
    for j=upInts
        [layerEqs,prob,nLayers]=numSubsIntUp(j,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,...
                               layerEqs,nLayers,...
                               Kx,b2,c1,c2,c5,...
                               termed,parents,inLayer);
    end

    for j=downInts
        [layerEqs,prob,nLayers]=numSubsIntDn(j,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,...
                               layerEqs,nLayers,...
                               Kx,b2,c1,c2,c5,...
                               termed,parents,inLayer);
    end
    
    %HERE add hookEqs to layerEqs
    layerEqs=cat(1,layerEqs,hookEqs);
    
    for j=prob.targ'
        if ismember(j,prob.tops)
            %take fclose and close it in all layers
            if any(iLinkClose==j)
            	layerEqs=numCloseTops(j,layerEqs,closeEqs,iLinkClose,nLayers);
            else
                layerEqs=numSubsIntUp(j,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,...
                                    layerEqs,nLayers,Kx,b2,c1,c2,c5,...
                                termed,parents,inLayer);
            end
        else
            layerEqs=numConnUD(j,layerEqs,closeEqs,iLinkClose,extraEqs,iLinkExtra,nLayers,prob,termed,Kx,parents);
            if termed(ismember(prob.iLinks,j))
                layerEqs=numCloseTops(j,layerEqs,closeEqs,iLinkClose,nLayers);
            end
        end
    end

    for i=1:nLayers
        topGradI=discoverIndices(layerEqs(i).vars,'G1');
        termGrad=ismember(topGradI,iLinkClose);
        if any(termGrad)
            for j=find(termGrad)'
                iClEq=iLinkClose==topGradI(j);
                layerEqs(i)=subsFor(layerEqs(i),closeEqs(iClEq).depvar,...
                    closeEqs(iClEq).vars,closeEqs(iClEq).coefs);
                layerEqs(i)=sumVars(layerEqs(i));
            end
        end
    end
    
    %need to also make extraEquations for targs that are on single strand
    %from top to bot, if bot is meant to be closed by probDef DOF
    %accounting...
    
    nDOF=sum(startsWith(cat(1,{layerEqs(:).depvar}),'psiXBar'))-1;
    allVars=unique(cat(2,layerEqs(:).vars));
    nVars=sum(startsWith(allVars,'psi1') | startsWith(allVars,'G1'));
    nH=sum(strcmp('H',cat(1,{layerEqs(:).kLayer})));
    n0=sum(strcmp('0',cat(1,{layerEqs(:).depvar})));
    if nDOF+nH+n0<nVars
        nExtraEqs=nVars-nDOF;
        %check that prob.targ is on unforked path from top to bot
        hasJunc=false(size(prob.targ));
        for i=1:size(prob.targ,1)
            hasJunc(i)=checkJuncs(prob.targ(i),prob.tops,prob.bots,parents);
        end
        if any(~hasJunc)
            canClose=prob.targ(~hasJunc);
            canClose=canClose(1:nExtraEqs);  %will error if there are other missing equations...
            nExtraClose=size(canClose,1);
            for i=1:nExtraClose
                j=distalTipSrch(canClose(i),parents);
                nTerms=size(j,1);
                [locCloseEqs,jLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
                [extraCloseEq,~]=numCloseToSeg(parents(canClose(i)),locCloseEqs,jLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
                for j=1:nLayers
                    if ismember(extraCloseEq.depvar,layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),extraCloseEq.depvar,...
                            extraCloseEq.vars,extraCloseEq.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                end
            end
        end
    end
end