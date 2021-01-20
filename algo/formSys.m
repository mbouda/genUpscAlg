function badVar=formSys(elimVars,eqs,fullSet,iEq)

    %currently used only to identify which variables are 'bad' and should
    %be taken out pre-solution
    %ideally, could expand this function to establish full order and method
    %of eliminating variables

    targEq=eqs(fullSet==iEq);
    remEqs=eqs(fullSet~=iEq);
    nRem=size(remEqs,1);
    
    allVars=unique(cat(2,eqs(:).vars));
    isElim=ismember(allVars,elimVars);
    allVars=cat(2,allVars(isElim),allVars(~isElim));
    
    inTarg=ismember(allVars,targEq.vars);
    nVars=size(allVars,2);
    
    %should have domainVars input into here, to use it to set up the set of
    %what should be in there and what should not... ?
    inRem=false(nRem,nVars);
    for i=1:nRem
        inRem(i,:)=ismember(allVars,remEqs(i).vars);
    end

    inEqs=cat(1,inRem,inTarg);
    
    badVar=allVars((sum(inEqs)==1 & ~inEqs(end,:) & isElim));
    
end