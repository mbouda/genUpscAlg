function eqs=solveSysFor(iEq,eqSet,fullSet,eqs,elimVars,nVars)

    eqList=cat(1,setdiff(eqSet,iEq),iEq);
    [~,eqInd]=ismember(eqList,fullSet);

    for j=1:nVars
        eqs(eqInd(j))=numIsolate(eqs(eqInd(j)),elimVars{j});
        for k=eqInd(j+1:end)'
            eqs(k)=subsFor(eqs(k),eqs(eqInd(j)).depvar,...
                            eqs(eqInd(j)).vars,...
                            eqs(eqInd(j)).coefs);
            eqs(k)=sumVars(eqs(k));
        end
    end

end