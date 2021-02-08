function eqSols=targEqSol(testCase,testSol)
keyboard

    Kr=pi*(2*testCase.params.r-testCase.params.b).*testCase.params.kr./testCase.params.b; 
    [~,CP1,CP2,CP3]=popJCPP(testCase.parents,...
                            Kr,testCase.params.Kx,testCase.L,testCase.nL);
    invCP2=inv(CP2);
    [~,CM1,CM2,CM3]=popJCMP(testCase.parents,...
                            Kr,testCase.params.Kx,testCase.L,testCase.nL);
    invCM2=inv(CM2);
    
    KrS=testCase.params.KrS;
    inLayer=testCase.inLayer;
        
    psi0=zeros(testCase.nL,testSol.nBC,testSol.nC);
    psiX=zeros(testCase.nL,testSol.nBC,testSol.nC);
    Qr=zeros(testCase.nL,testSol.nBC,testSol.nC);
    

    for j=1:testSol.nC
        psiC=testSol.psiBC(j);
        parfor i=1:testSol.nBC
            psiSL=testSol.psiBCS(:,i);
            psiS=psiSL(inLayer);
            psi0(:,i,j)=-invCP2*(CP3*psiS+psiC*CP1);
            psiX(:,i,j)=-invCM2*(CM3*psiS+psiC*CM1);
            Qr(:,i,j)=-KrS.*(psiX(:,i,j)-psiS);
        end
    end
  
    psi1=zeros(testCase.nL,testSol.nBC,testSol.nC);
    psi1(1,:)=psiC;
    psi1(2:end,:)=psi0(testCase.parents(2:end),:);
    
    G1=zeros(testCase.nL,testSol.nBC,testSol.nC);
    for j=1:testSol.nC
        parfor i=1:testSol.nBC
            psiSL=testSol.psiBCS(:,i);
            psiS=psiSL(inLayer);
            G1(:,i,j)=(psiX(:,i,j)-testCase.params.c1.*psi1(:,i,j)-(1-testCase.params.c1).*psiS)...
                ./testCase.params.c2;
        end
    end
    
    for i=1:size(testCase.prob,1)
        
        eqs=testCase.prob(i).eqs;
        nEqs=size(eqs,1);
        eqSols=struct('vars',cell(nEqs,1),'sols',cell(nEqs,1));

        for j=1:nEqs

            vals=zeros(testSol.nBC,size(eqs(j).vars,2));
            isPsiX=startsWith(eqs(j).vars,'psiXBar');
            isPsiL=startsWith(eqs(j).vars,'psiL');
            isBC=endsWith(eqs(j).vars,'C');
            isPsi1=startsWith(eqs(j).vars,'psi1');
            isG1=startsWith(eqs(j).vars,'G1');

            jL=discoverIndices(eqs(j).vars(isPsiL),'psiL');
            jX=discoverIndices(eqs(j).vars(isPsiX),'psiXBar');
            jp1=discoverIndices(eqs(j).vars(isPsi1),'psi1');
            jg1=discoverIndices(eqs(j).vars(isG1),'G1');
            jC=discoverIndices(eqs(j).vars(isBC),'psiC');  %currently not working for G1...

            
            
            psiLj=testSol.psiBCS(jL,:)';
            psi1j=psi1(jp1,:,1)';
            G1j=G1(jg1,:,1)';
            
            vals(:,isPsiL)=psiLj;
            vals(:,isPsi1)=psi1j;
            vals(:,isG1)=G1j;
 
            %for main layerEqs:
            depJ=discoverIndices({eqs(j).depvar},'psiXBar');
            depVal=testSol.psiXL(depJ,:)';

            %forHookEqs:
            depJ=discoverIndices({eqs(j).depvar},'psi1');
            depVal=psi1(depJ,:,1)';
            
            %for zeroEqs:
            depVal=0;
            
            %for new extraEqs
            depJ=discoverIndices({eqs(j).depvar},'G1');
            depVal=G1(depJ,:,1)';
            
            predVal=sum(repmat(eqs(j).coefs,[testSol.nBC,1]).*vals,2);
            r=depVal-predVal;
            
            figure; histogram(r)
            figure; plot(depVal,predVal,'.')
            %these should then be compared to solutions made with coefs..
            
            pred2Val=sum(repmat(KrS(inLayer==depJ),[1 testSol.nBC]).*psiX(inLayer==depJ,:))'./sum(KrS(inLayer==depJ));
            
            figure; plot(predVal,pred2Val,'.')
            figure; plot(depVal,pred2Val,'.')
            
            %can just use testSol to calculate without re-solving?
                %no, because have psiX for layers, not segments
            
            testCase
            
            testSol
            
            
        end

    end

end