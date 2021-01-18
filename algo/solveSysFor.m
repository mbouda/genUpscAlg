function eqs=solveSysFor(iEq,eqSet,fullSet,eqs,elimVars,nVars)


    eqList=setdiff(eqSet,iEq);
    pres=false(nVars);
    for i=1:nVars
        for j=1:nVars
            pres(i,j)=ismember(elimVars{i},eqs(fullSet==eqList(j)).vars);
        end
    end
    
    badVar=~any(pres,2);
    if any(badVar)
        remVars=unique(cat(2,eqs(ismember(fullSet,eqList)).vars));
        iL=discoverIndices(remVars,'psiL');
        iBV=find(badVar);
        for i=iBV'
            %exchange for psiL from the others
            [~,jGV]=max(abs(iL-iEq));
            elimVars{i}=sprintf('psiL%d',iL(jGV));
            iL(jGV)=[];
            for j=1:nVars
                pres(i,j)=ismember(elimVars{i},eqs(fullSet==eqList(j)).vars);
            end
        end
    end
    
    nEqVar=sum(pres,2);
    [~,iVar]=sort(nEqVar);
    eqI=zeros(nVars,1);
        
    for i=iVar(1:end-1)'
        inEqs=find(pres(i,:));
        elEqs=sum(pres(:,pres(i,:)),1)~=1;
        eqI(i)=inEqs(find(elEqs,1));
        pres(:,eqI(i))=false;
    end
    eqI(iVar(end))=find(pres(iVar(end),:)); %should be 1 left, only
    
    jEq=fullSet==iEq;
    for i=1:nVars
        EQi=fullSet==eqList(eqI(i));
        eqs(EQi)=numIsolate(eqs(EQi),elimVars{i});
        for j=i+1:nVars
            eqJ=fullSet==eqList(eqI(j));
            if ismember(elimVars{i},eqs(eqJ).vars)
                eqs(eqJ)=subsFor(eqs(eqJ),eqs(EQi).depvar,...
                            eqs(EQi).vars,eqs(EQi).coefs);
                eqs(eqJ)=sumVars(eqs(eqJ));
            end
        end
        if ismember(elimVars{i},eqs(jEq).vars)
            eqs(jEq)=subsFor(eqs(jEq),eqs(EQi).depvar,...
                        eqs(EQi).vars,eqs(EQi).coefs);
            eqs(jEq)=sumVars(eqs(jEq));
        end
    end
end
