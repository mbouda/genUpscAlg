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
    nBL=prob.iLinks(nonBase(prob.iLinks));
    
    exts=cat(1,prob.iLinks(~nonBase(prob.iLinks)),...
        nBL(~ismember(inLayer(parents(nBL)),prob.kLayers)),...
        parents(dtr(~ismember(inLayer(ismember(parents,prob.iLinks)),prob.kLayers))));

    prob.ints=setdiff(prob.iLinks,cat(1,exts,prob.terms));
    dtr=find(ismember(parents,exts));
    par=parents(dtr);
    prob.bots=unique(par(ismember(inLayer(par),prob.kLayers) & ~ismember(inLayer(dtr),prob.kLayers)));
    tf=false(size(parents));
    tf(nonBase)=~ismember(inLayer(parents(nonBase)),prob.kLayers);
    prob.tops=exts(~nonBase(exts) | (ismember(inLayer(exts),prob.kLayers) & tf(exts)));
    
    sectors=pickSectors(prob.tops,prob.iLinks,parents);
    nSect=size(sectors,1);
    
    if iLayer==1
        prob.targ=find(parents==0);
    else
        notBase=parents>0;
        targs=false(size(parents));
        targs(notBase)=inLayer(notBase)==iLayer & inLayer(parents(notBase))~=iLayer;
        
        %here, test whether they descend from another targ or a top
            %if targ, make false
        iTarg=find(targs);
        isDesc=srchTargPars(iTarg,prob.tops,parents);
        targs(iTarg)=~isDesc;
        nTarg=floor((bot-top+all(parents(prob.tops)==0)+all(ismember(prob.bots,prob.terms))+...
            (size(prob.tops,1)-1)-(size(prob.tops,1)-size(unique(parents(prob.tops)),1)))/2);
        
        %take floor on assumption that if get 1.5, can't fully close that extra link
        %if all tops are bases, then can eliminate extra dof with collar condition
        %if all bots are terms, then eliminate extra dof b/c netowrk is closed
            %note: in this case, bots is actually empty, all([]) -> true
        %extra tops mean multiple separate sub-networks; need at least one
        %target per sub-network, right?
            
        prob.targ=find(targs,nTarg,'first');
        nTargs=sum(targs);
        if nTargs<nTarg
            %can add difference from different sectors
            targSect=zeros(nTargs,1);
            for i=1:nTargs
                j=1;
                while ~ismember(prob.targ(i),sectors{j})
                    j=j+1;
                end
                targSect(i)=j;
            end
            
            remnSect=setdiff((1:nSect)',targSect);
            if any(remnSect)
                nAdd=nTarg-nTargs;
                for i=1:nAdd
                    canBeTarg=inLayer(sectors{remnSect(i)})==iLayer & ...
                        inLayer(parents((sectors{remnSect(i)})))~=iLayer;
                    if any(canBeTarg)
                        prob.targ=cat(1,prob.targ,sectors{remnSect(i)}(find(canBeTarg,1)));
                    elseif any(ismember(sectors{remnSect(i)},prob.tops))
                        prob.targ=cat(1,prob.targ,sectors{remnSect(i)}(ismember(sectors{remnSect(i)},prob.tops)));
                    else
                        crosses=inLayer(sectors{remnSect(i)})~=inLayer(parents((sectors{remnSect(i)})));
                        prob.targ=cat(1,prob.targ,sectors{remnSect(i)}(find(crosses,1)));
                    end
                end
            end
        end
    end
    
    
end