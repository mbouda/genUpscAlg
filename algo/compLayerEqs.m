function layerEqs=compLayerEqs(iLayer,layerEqs,prob,b2,c1,c2,c5,Kx,inLayer,parents,nLayers)

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
            
    % ints in order...
    downInts=intersect(prob.ints(inLayer(prob.ints)<iLayer),prob.iLinks(~termed))';
    upInts=flipud(intersect(prob.ints(inLayer(prob.ints)>=iLayer),prob.iLinks(~termed)))';
    
    for j=1:nOJ
        nLegsTarg=numel(extraEqs(j).targTrack);
        if nLegsTarg==1
            downExtra=ismember(downInts,iLinkExtra(j));
            downExtra=runDownExtras(iLinkExtra(j),downExtra,downInts,parents);
            topExtra=ismember(prob.tops,iLinkExtra(j));
            downExtra=runDownExtras(iLinkExtra(j),downExtra,downInts,parents);
            upInts=cat(2,upInts,prob.tops(topExtra));
            upInts=cat(2,upInts,downInts(downExtra));
            downInts=downInts(~downExtra);
            sib=setdiff(cat(1,extraEqs(j).helperEqs(:).iLink),extraEqs(j).iLink);
            if ~ismember(sib,downInts)
                downInts=cat(2,downInts,sib);
            end
        elseif nLegsTarg==2
            keyboard
            %will need to be both up AND down!
        end %if ==0, no need to execute anything
    end
    downInts=cat(2,downInts,topTermSibs);
    
    upInts=sort(upInts,'descend');
    downInts=sort(downInts,'ascend');

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
    
    %link across top interface of layer i
%         layerSegs=find(inLayer==iLayer);
%         layerSegs=layerSegs(inLayer(parents(layerSegs))~=iLayer);
        %the above two lines ended up too permissive, b/cof the closeTops
        %on termed below
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

end