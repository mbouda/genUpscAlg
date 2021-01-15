function outs=defLayerEqs(nLay,kLayers,inLayers,KrS)

    outs=struct('kLayer',cell(nLay,1),...
        'coefs',cell(nLay,1),'vars',cell(nLay,1),...
        'depvar',cell(nLay,1));
    for i=1:nLay
        outs(i).kLayer=kLayers(i);
        iSeg=inLayers==kLayers(i);
        outs(i).coefs=KrS(iSeg)'/sum(KrS(iSeg));
        nSeg=sum(iSeg);
        jSeg=find(iSeg);
        varNames=cell(1,nSeg);
        for j=1:nSeg
            varNames{j}=sprintf('psiX%d',jSeg(j));
        end
        outs(i).vars=varNames;
        outs(i).depvar=sprintf('psiXBar%d',kLayers(i));
    end

end