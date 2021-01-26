function downLinks=idDownLinks(targs,tops,parents)

    nTargs=size(targs,1);
    lists=cell(nTargs,1);
    for i=1:nTargs
        linkList=[];
        if ~ismember(targs(i),tops)
            l=parents(targs(i));
            while ~ismember(l,tops) % && l>0 %should not be necessary
                linkList=cat(1,linkList,l);
                l=parents(l);
            end
            lists{i}=linkList;
        end
    end
    downLinks=cat(1,lists{:});    
end