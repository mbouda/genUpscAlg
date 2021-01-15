function check=checkSol(testCase,fullSol,tol,resTol)

    for i=1:testCase.nDomLayers
        
        solLayers=discoverIndices(testCase.sol(i).vars,'psiL');
        othLayers=discoverIndices(testCase.sol(i).vars,'psiXBar');
        
        isBdry=ismember('psiC',testCase.sol(i).vars);
        
        warning off   %this can likely be improved, by (a) specifying which, (b) storing warning in check structure...
        c=[fullSol.psiBCS(solLayers,:)' fullSol.psiXL(othLayers,:)' repmat(fullSol.psiBC,[size(fullSol.psiXL,2) isBdry])]\fullSol.psiXL(i,:)';
        warning on
        C=sort([c testCase.sol(i).coefs']);
        
        check(i).hard=C(:,1)==C(:,2);
        
        check(i).soft=abs((C(:,1)-C(:,2))./C(:,1))<tol | abs(C(:,1))<tol;
        check(i).allSoft=all(check(i).soft);
        
        predSol=zeros(fullSol.nBC,1);
        for j=1:size(solLayers,1)
        	predSol=predSol+fullSol.psiBCS(solLayers(j),:)'*testCase.sol(i).coefs(ismember(testCase.sol(i).vars,sprintf('psiL%d',solLayers(j))));
        end
        for j=1:size(othLayers,1)
        	predSol=predSol+fullSol.psiXL(othLayers(j),:)'*testCase.sol(i).coefs(ismember(testCase.sol(i).vars,sprintf('psiXBar%d',othLayers(j))));
        end
        for j=1:sum(isBdry)
            predSol=predSol+fullSol.psiBC*testCase.sol(i).coefs(ismember(testCase.sol(i).vars,sprintf('psiC')));
        end

        r=(predSol-fullSol.psiXL(i,:)')./fullSol.psiXL(i,:)';
        check(i).maxRes=max(abs(r));
        check(i).resTest=check(i).maxRes<resTol;
    end
    
end

  %figure; histogram(r)
  %figure; plot(fullSol.psiXL(i,:)',predSol,'.')