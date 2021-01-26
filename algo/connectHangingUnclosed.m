function  hookEqs=connectHangingUnclosed(hangLinks,prob,parents,b2,c1,c2,c5,Kx,inLayer)

    nHL=size(hangLinks,1);
    for i=1:nHL
        j=distalTipSrch(hangLinks(i),parents);
        nTerms=size(j,1);
        [locCloseEqs,jLinkClose]=numCloseTerms(j,b2(j),c1(j),c2(j),inLayer(j),nTerms);
        [locCloseEqs,jLinkClose]=numCloseToSeg(parents(hangLinks(i)),locCloseEqs,jLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer);
    end
    
    hookEqs=connectHanging(hangLinks,locCloseEqs,jLinkClose,prob,parents,b2,c1,c2,c5,Kx,inLayer);
            
end