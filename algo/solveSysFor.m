function eqs=solveSysFor(iEq,eqSet,fullSet,eqs,elimVars,nVars)

%need to add provision for depvar==0
% including fact, that 0 is not in fullSet!!
    %vars not in iEq need not be considered (?) this may assume too much
    
    mLayers=size(eqs,1);
    nLayers=size(fullSet,1); %should bring in from outside...
    if mLayers>nLayers
        fullSet=cat(1,eqs(:).kLayer);
        eqSet=fullSet;
    end
    
    
    badVars=formSys(elimVars,eqs,fullSet,iEq);
    
%     relVars=intersect(elimVars,eqs(fullSet==iEq).vars);
%     badVars=setdiff(elimVars,relVars);
    
    nVars=size(elimVars,2);
    
    eqList=setdiff(eqSet,iEq);
    nEqs=size(eqList,1);
    pres=false(nVars,nEqs);
    for i=1:nVars
        for j=1:nEqs
            pres(i,j)=ismember(elimVars{i},eqs(fullSet==eqList(j)).vars);
        end
    end
    
    for i=1:size(badVars,2)
        isBad=ismember(elimVars,badVars{i})';
        remEq=find(pres(isBad,:),1);
        pres(:,remEq)=false;
        pres(isBad,:)=[];
        elimVars(isBad)=[];
        nVars=nVars-1;
    end
    %take out one equation per badVar...
    
    
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
    eqI(iVar(end))=find(pres(iVar(end),:),1);
    
    jEq=fullSet==iEq;
    for i=1:nVars
        EQi=fullSet==eqList(eqI(i));
        if eqs(EQi).kLayer==0
            eqs(EQi)=numIsolateInZero(eqs(EQi),elimVars{i});
        else
            eqs(EQi)=numIsolate(eqs(EQi),elimVars{i});
        end
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
