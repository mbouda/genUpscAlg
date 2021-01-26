function sec=growSec(l,sec,parents,iLinks)

    if ismember(l,iLinks)
        sec=cat(1,sec,l);
        dtrs=find(parents==l);
        for d=dtrs'
            sec=growSec(d,sec,parents,iLinks);
        end
    end

end