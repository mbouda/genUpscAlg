function downExtra=runDownExtras(p,downExtra,downInts,parents)

    dtrs=find(parents==p);
    nDtrs=size(dtrs,1);
    for i=1:nDtrs
        downDtr=ismember(downInts,dtrs(i));
        if any(downDtr)
            downExtra=runDownExtras(dtrs(i),downExtra,downInts,parents);
            downExtra = downExtra | downDtr;
        end
    end
end