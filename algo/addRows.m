function cM=addRows()

keyboard
    %need to adapt this code to its new position...

    %% below this, meant to fix systems with subsystems of fewer variables than equations
        %but seems to break more cases than helps
        %likely should integrate it into formSys, and supply necessary
        %inputs to that...
        %there could check if can also eliminate some equations, all at once
        %and ideally arrive at the best-shaped system...
    %put this stuff below into standalone function...
    allVars=unique(cat(2,layerEqs(:).vars));
    N=size(allVars,2);
    
    isPsiX=startsWith(allVars,'psiXBar');
    isPsiL=startsWith(allVars,'psiL');
    isBC=endsWith(allVars,'C');
    isNot=~(isPsiX | isPsiL | isBC);
    
    depVars=cat(2,{layerEqs(:).depvar});
    remEq=ismember(depVars,sprintf('psiXBar%d',iLayer)) | ~(startsWith(depVars,'psiX') | startsWith(depVars,'psiL') | endsWith(depVars,'C'));
    %remove also the psiL equations?
    
    
    depVars=unique(depVars);
    
    toElim=setdiff(allVars(isNot),depVars);
    %should also check if they're in the target equation,. no? don't bother
    %if not(?)
    
    nElim=size(toElim,2);
    nEqs=size(layerEqs,1);
    varInEq=false(nEqs,nElim);
    for i=1:nElim
        for j=1:nEqs
            varInEq(j,i)=ismember(toElim{i},layerEqs(j).vars);
        end
    end
    varInEq(remEq',:)=[];

    badCombs=cell(0,1);
    for i=1:nElim
        if sum(varInEq(:,i))<nElim
            oth=setdiff(1:nElim,i);
            for j=1:nElim-1
                set=false(1,size(oth,2));
                set(1:j)=true;
                combs=unique(perms(set),'rows');
                ncombs=size(combs,1);
                for k=1:ncombs
                    if sum(varInEq(:,i) | varInEq(:,oth(combs(k,:))))<j+1
                        newComb=sort([i oth(combs(k,:))]);
                        nC=numel(badCombs);
                        if nC>0
                            %only add if separate from all existing ones
                            
                            isNew=true;
                            for l=1:nC
                                isNew=isNew & ~all(ismember(newComb,badCombs{l}));
                            end
                            if isNew
                                badCombs{end+1,1}=[i oth(combs(k,:))];
                            end
                        else
                            badCombs{end+1,1}=newComb;
                        end
                    end
                end
            end
        end
    end
    
    nC=size(badCombs,1);
    if nC>1
        for i=1:nC
            oC=setdiff(1:nC,i);
            isCont=false(nC-1,1);
            for j=oC
                isCont(j)=all(ismember(badCombs{i},badCombs{j}));
            end
            if any(isCont)
                badCombs{i}=[];
            end
        end
        remComb=false(nC,1);
        for i=1:nC
            remComb(i)=numel(badCombs{i})==0;
        end
        badCombs(remComb)=[];
        nC=size(badCombs,1);
    end
    
    hangVars=cell(1,0);
    for i=1:nC
        hangVars=cat(2,hangVars,toElim(badCombs{i}(1:end-1)));
    end
    psiHang=discoverIndices(hangVars,'psi1');
    gHang=discoverIndices(hangVars,'G1');
    if any(gHang)
        warning('Hanging G1','hanG');
        %keyboard
        %if legitimately here, would need another king of hook equation,
        %since for now, they work for closed links
    end
    if any(psiHang)
        hookEqs=connectHangingTarg(psiHang,closeEqs,iLinkClose,prob,parents,b2,c1,c2,c5,Kx,inLayer);
        hookEqs=rmfield(hookEqs,'iLink');
        [hookEqs(:).kLayer]=deal('H');
        layerEqs=cat(1,layerEqs,hookEqs);
    end
    
    

end