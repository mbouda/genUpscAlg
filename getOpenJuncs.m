function ojps=getOpenJuncs(iLinks,parents,termed)

    pars=ismember(iLinks,parents);
    isOpen=false(size(pars));
    for j=iLinks(pars)'
        dtrs=parents==j;
        if sum(dtrs)>1
            isOpen(iLinks==j)=sum(~termed(ismember(iLinks,find(dtrs))))>1;
        end
    end
    ojps=iLinks(isOpen);

    %find sibs that are both untermed
    nonDomPar=setdiff(parents(iLinks),iLinks);
    isOpen=false(size(nonDomPar));
    for j=nonDomPar'
        dtrs=parents==j;
        if sum(dtrs)>1
            isOpen(nonDomPar==j)=sum(~termed(ismember(iLinks,find(dtrs))))>1;
        end
    end
    ojps=cat(1,ojps,nonDomPar(isOpen));

end