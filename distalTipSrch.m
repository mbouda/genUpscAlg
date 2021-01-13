function distTips=distalTipSrch(j,parents)

    isTip=false(size(parents));
    dtrs=parents==j;
    isTip=findTips(dtrs,parents,isTip);
    distTips=find(isTip);
    
end


