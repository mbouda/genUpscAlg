function [inLay,sameLay,crossLay,nLay]=linkLayersC(cz,zLims,nL,nLay)

    zEnds=[cz(:,2) sum(cz,2)];

    inLay=zeros(nL,2);
    for i=1:2
        lt=repmat(zEnds(:,i),[1 nLay+1])<repmat(zLims',[nL 1]);
        lt(:,1)=lt(:,1)|zEnds(:,i)==zLims(1); %replace first row w/ >=, to place u=u_0 on first interval
        [inLay(:,i),~]=find((lt(:,1:nLay) & ~lt(:,2:nLay+1))'); %fails here if any segments are above 0
        %[~,jj]=sort(jj);
        %inLay(:,i)=inLay(jj,i);
    end
    
    sameLay=inLay(:,1)==inLay(:,2);
    crossLay=~sameLay;
    
end