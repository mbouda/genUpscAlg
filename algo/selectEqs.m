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
        nEC=zeros(nC,1);
        for j=1:nC
            nEC(j)=sum(cands(j,:));
            nVC(j)=sum(elimCoef | cands(j,:));
        end
        
        hasLeastV=nVC==min(nVC);
        iL=find(hasLeastV);
        if sum(hasLeastV)>1
            hasMostE=nEC(hasLeastV)==max(nEC(hasLeastV));
            if sum(hasMostE)>1
                iE=find(hasMostE);
                
                xCDom=sum(domCoef(hasMostE,:),2);
                [~,ice]=min(xCDom);
                icd=iE(ice);
            else
                icd=find(hasMostE);
            end
        else
            icd=1;
        end
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
    
    %HERE: if can get rid of further variables by adding equations, without
    %adding isNot variables, then should?
    
    cM=cM(nT-nE+1:nT,:);
    
    cut=isNot & ~any(cM);
    cM(:,cut)=[];

end