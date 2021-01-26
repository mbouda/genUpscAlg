function isDesc=srchTargPars(targs,tops,parents)

    nTargs=size(targs,1);
    isDesc=false(nTargs,1);
    for i=1:nTargs
        l=parents(targs(i));
        while ~ismember(l,targs) && ~ismember(l,tops)
            l=parents(l);
        end
        isDesc(i)=ismember(l,targs);
    end

end