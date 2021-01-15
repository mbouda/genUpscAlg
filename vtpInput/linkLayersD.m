function inLayer=linkLayersD(zEnds,zLims,nL,nLay)

    [endIsLim,interfI]=ismember(zEnds,zLims);
    
    inLayer=zeros(nL,1);
    for i=1:nL
        intersectLayers=cell(1,2);
        for j=1:2
            if endIsLim(i,j)
                intersectLayers{j}=[interfI(i,j)-1 interfI(i,j)];
            else
                lt=repmat(zEnds(i,j),[1 nLay+1])<repmat(zLims',[nL 1]);
                [intersectLayers{j},~]=find((lt(:,1:nLay) & ~lt(:,2:nLay+1))');
            end
        end
        inLayer(i)=intersect(intersectLayers{1},intersectLayers{2});
    end
        
end