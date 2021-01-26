function sectors=pickSectors(tops,iLinks,parents)

    nTops=size(tops,1);
    sectors=cell(nTops,1);
    for i=1:nTops
        sectors{i}=growSec(tops(i),[],parents,iLinks);
    end
    
    if any(setdiff(iLinks,cat(1,sectors{:})))
        warning('not all links sectored out','extrL');
    end
    
end