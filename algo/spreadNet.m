function subNet=spreadNet(l,subNet,parents,isSub)

    p=parents(l);
    d=find(parents==l);
    s=setdiff(find(parents==p),l);
    
    if any(s) %necessary b/c setdiff screws up dimension;
        J=cat(1,p(p>0),s,d);
    else
        J=cat(1,p(p>0),d);
    end
    jSub=isSub(J);
    K=setdiff(J(jSub),subNet);
    
    if any(K)
        subNet=cat(1,subNet,K);
        for k=K'
            subNet=spreadNet(k,subNet,parents,isSub);
        end
    end
end
