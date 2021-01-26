function hasJunc=checkJuncs(targ,tops,bots,parents)

    hasJunc=false;
    l=targ;
    while ~hasJunc && ~ismember(l,tops)
        l=parents(l);
        hasJunc=sum(parents==l)>1;
    end
    
    l=targ;
    while ~hasJunc && ~ismember(l,bots)
        hasJunc=sum(parents==l)>1;
        l=find(parents==l);
    end
end