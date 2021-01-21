function [cM,cut]=selectEqs(cM,isNot,nT)

    nDom=sum(~isNot);

    hasCoef=cM~=0;
    domCoef=hasCoef(:,~isNot).*repmat(nDom:-1:1,[nT 1]);

    nE=1;
    elimCoef=hasCoef(end,isNot);
    nV=sum(elimCoef);
    
    while nV>=nE
        nC=nT-nE;
        cands=hasCoef(1:nC,isNot);
        nVC=zeros(nC,1);
        for j=1:nC
            nVC(j)=sum(elimCoef | cands(j,:));
        end
        hasLeastE=nVC==min(nVC);
        if sum(hasLeastE)>1
            xCDom=sum(domCoef(hasLeastE,:),2);
            [~,icd]=min(xCDom);
        else
            icd=1;
        end
        iL=find(hasLeastE);
        iEq=iL(icd);
        
        in=cM(iEq,:);
        out=cM(nC,:);
        
        cM(nC,:)=in;
        cM(iEq,:)=out;
        
        hasCoef=cM~=0;
        domCoef=hasCoef(:,~isNot).*repmat(nDom:-1:1,[nT 1]);
        
        nV=nVC(iEq);
        nE=nE+1;
    end
    
    cM=cM(nT-nE+1:nT,:);
    
    cut=isNot & ~any(cM);
    cM(:,cut)=[];

end