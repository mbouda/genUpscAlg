function  sol=matSol(sys,psiBC,psiBCS,nBC,nC,nL)

    nS=size(sys,2);

    sol=struct('psiXL',cell(nS,1));
    
    for iS=1:nS
        
        mySys=sys(iS);
        [~,CB1,CB2,CB3,~,~]=popJCMP(mySys.parents,mySys.Kr,mySys.Kx,mySys.L,mySys.nL);
        invCB2=inv(CB2);
        
        KrS=mySys.Kr.*mySys.L; 
        inLayer=mySys.inLayer;
        
        psiXS=zeros(mySys.nL,nBC,nC);
        Qr=zeros(mySys.nL,nBC,nC);

        for j=1:nC
            psiC=psiBC(j);
            parfor i=1:nBC
                psiSL=psiBCS(:,i);
                psiS=psiSL(inLayer);
                psiXS(:,i,j)=-invCB2*(CB3*psiS+psiC*CB1);
                Qr(:,i,j)=-KrS.*(psiXS(:,i,j)-psiS);
            end
        end
        
        D=repmat(KrS,[1 nBC nC]);
        N=psiXS.*D;
        sol(iS).psiXL=zeros(nL,nBC,nC);
        sol(iS).QR=zeros(nL,nBC,nC);
        for iL=1:nL
            isInLayer=inLayer==iL;
            sol(iS).psiXL(iL,:,:)=sum(N(isInLayer,:,:),1)./sum(D(isInLayer,:,:),1);
            sol(iS).QR(iL,:,:)=sum(Qr(isInLayer,:,:),1);
        end
        
    end

end