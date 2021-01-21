function prob=numProbDef(iLayer,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer)

    top=iLayer;
    bot=iLayer;
    low=max(inLayer);  %bring in from outside!
    
    nDOF=subDomDOFs(top,bot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
    
    while (bot-top)<nDOF
        provBot=bot+1;
        if provBot>low
            provBot=low;
        end
        provTop=top-1;
        if provTop<1
            provTop=1;
        end
        
        nDOFDn=subDomDOFs(top,provBot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
        nDOFUp=subDomDOFs(provTop,bot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
        nDOFUD=subDomDOFs(provTop,provBot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
        
        
        %for current case, not optimal: includes 2 b/c closest to 3, even
        %though can do from layers 3--5
        if nDOFDn<=(provBot-top) && top+1<=iLayer
            %test reduction by 1 from top:
            nDOFt=subDomDOFs(top+1,provBot,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
            if nDOFt<=(provBot-(top+1))
                nDOFDn=nDOFt;
                top=top+1;
            end
            
            bot=provBot;
            nDOF=nDOFDn;
        elseif nDOFUp<=(bot-provTop) && bot-1>=iLayer
            %test reduction by 1 from bot:
            nDOFt=subDomDOFs(provTop,bot-1,acroLrs,nxtLr,basiLrs,prvLr,parents,inLayer);
            if nDOFt<=((bot-1)-provTop)
                nDOFUp=nDOFt;
                bot=bot-1;
            end
            
            top=provTop;
            nDOF=nDOFUp;
        else
            top=provTop;
            bot=provBot;
            nDOF=nDOFUD;
        end
        
    end
    
    prob.iLinks=find(inLayer>=top & inLayer<=bot);
    prob.kLayers=(top:bot)';
    
    isPar=ismember(prob.iLinks,parents); 
    prob.terms=prob.iLinks(~isPar);  %identifies terminals

    dtr=find(ismember(parents,prob.iLinks));
    
    nonBase=parents~=0;
    
    
    exts=cat(1,prob.iLinks(~nonBase(prob.iLinks)),...
        prob.iLinks(~ismember(inLayer(parents(prob.iLinks(nonBase(prob.iLinks)))),prob.kLayers)),...
        parents(dtr(~ismember(inLayer(ismember(parents,prob.iLinks)),prob.kLayers))));

    prob.ints=setdiff(prob.iLinks,cat(1,exts,prob.terms));
    %separate exts into tops and bots (for now, assumes unidirectional growth:)
    dtr=find(ismember(parents,exts));
    par=parents(dtr);
    prob.bots=unique(par(inLayer(par)==bot & inLayer(dtr)>bot));
    
    tf=false(size(parents));
    tf(nonBase)=inLayer(parents(nonBase))<top;
    prob.tops=exts(~nonBase(exts) | (inLayer(exts)==top & tf(exts)));

    if iLayer==1
        prob.targ=find(parents==0);
    else
        notBase=parents>0;
        targs=false(size(parents));
        targs(notBase)=inLayer(notBase)==iLayer & inLayer(parents(notBase))~=iLayer;
        
        
        nTarg=floor((bot-top+all(parents(prob.tops)==0)+all(ismember(prob.bots,prob.terms))+...
            (size(prob.tops,1)-1)-(size(prob.tops,1)-size(unique(parents(prob.tops)),1)))/2);
        %take floor on assumption that if get 1.5, can't fully close that extra link
        %if all tops are bases, then can eliminate extra dof with collar condition
        %if all bots are terms, then eliminate extra dof b/c netowrk is closed
            %note: in this case, bots is actually empty, all([]) -> true
        %extra tops mean multiple separate sub-networks; need at least one
        %target per sub-network, right?
        
            
        prob.targ=find(targs,nTarg,'first');
    end
    
    
end