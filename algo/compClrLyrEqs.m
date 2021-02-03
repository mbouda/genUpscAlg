function layerEqs=compClrLyrEqs(iLayer,collarCond,layerEqs,prob,b2,c1,c2,c5,Kx,inLayer,parents,nLayers)

    layerSegs=find(inLayer==iLayer);
    baseSegs=layerSegs(parents(layerSegs)==0);
    layerSegs=union(layerSegs(inLayer(parents(setdiff(layerSegs,baseSegs)))~=iLayer),baseSegs);

    topTerms=intersect(prob.terms,prob.tops);  %may always be empty here... eliminate for efficiency?
    topTermSibs=[];
    if any(topTerms)
        prob.tops(prob.tops==topTerms)=[];
        for j=topTerms'
            topTermSibs=cat(2,topTermSibs,setdiff(find(parents==parents(j)),j));
        end
    end
    
    j=setdiff(prob.tops,baseSegs);
    topLK=ismember(prob.kLayers,inLayer(j));
    [layerEqs(topLK),lidEqs]=subsCollar(j,collarCond,layerEqs(topLK),c1(j),c2(j),c5(j),b2(j),prob.kLayers(topLK),inLayer(j),numel(j));
    iLidLinks=cat(1,lidEqs(:).iLink);
    for j=find(topLK)'
            layerEqs(j)=sumVars(layerEqs(j));
    end
    
    j=intersect(prob.tops,baseSegs);
    topLK=ismember(prob.kLayers,inLayer(j));
    layerEqs(topLK)=subsBase(j,collarCond,layerEqs(topLK),c1(j),c2(j),prob.kLayers(topLK),inLayer(j),numel(j));
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
%     downInts=intersect(prob.ints(inLayer(prob.ints)<iLayer),prob.iLinks(~termed))';
%     upInts=flipud(intersect(prob.ints(inLayer(prob.ints)>=iLayer),prob.iLinks(~termed)))';

    allTops=union(prob.tops,topTerms);
    if size(allTops,2)>1
        allTops=allTops';
    end
    downLinks=idDownLinks(prob.targ,allTops',parents);
    downInts=intersect(downLinks,prob.ints)';
    upInts=setdiff(setdiff(prob.ints,downInts),prob.iLinks(termed))';
    
    for j=1:nOJ
        nLegsTarg=numel(extraEqs(j).targTrack);
        if nLegsTarg==1
%             downExtra=ismember(downInts,iLinkExtra(j));
%             downExtra=runDownExtras(iLinkExtra(j),downExtra,downInts,parents);
%             upInts=cat(2,upInts,downInts(downExtra));
%             downInts=downInts(~downExtra);
            
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
            
            %need to implement here the code from compLayerEqs?
                %(1) should make the definitive version a funciton and use
                %it both places
                %(2) if have collar, the situation of out-of-domain parents
                %shouldn't happen...
            
            %will any ever need passing BOTH up and down?
            %if so, in what order?
            
%         elseif nLegsTarg==2
%             keyboard
            %this option appears not no need any special code...
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
        jLidEq=iLidLinks==parents(j);
        if any(jLidEq)
            lidEqs(jLidEq)=numPassLidDn(j,prob,closeEqs,iLinkClose,extraEqs,iLinkExtra,lidEqs(jLidEq),Kx,b2,c1,c2,c5,termed,parents,inLayer);
            iLidLinks=cat(1,lidEqs(:).iLink);
        end
    end

    %after this, some of the extra Layer equations may remain unresolved

    %link across top interface of layer i

    for j=prob.targ'
        if ismember(j,baseSegs)
            [layerEqs,prob,nLayers]=numSubsBase(j,prob,collarCond,closeEqs,iLinkClose,extraEqs,iLinkExtra,layerEqs,nLayers,Kx,b2,c1,c2,c5,termed,parents,inLayer);
        else
            [layerEqs,prob,nLayers]=numConnUDClr(j,layerEqs,closeEqs,iLinkClose,lidEqs,extraEqs,iLinkExtra,collarCond,nLayers,prob,termed,Kx,parents);
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
end