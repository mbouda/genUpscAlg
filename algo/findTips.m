function isTip=findTips(dtrs,parents,isTip)
    
    J=find(dtrs)';
    for j=J
        dtrs=parents==j;
        if any(dtrs)
            isTip=findTips(dtrs,parents,isTip);
        else
            isTip(j)=true;
        end
    end

end