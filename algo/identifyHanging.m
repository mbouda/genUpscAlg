function hangLinks=identifyHanging(links,prob,parents)

    nLinks=size(links,1);
    lands=false(nLinks,1);
    for i=1:nLinks
        lands(i)=srchTarg(links(i),parents,prob.targ,cat(1,prob.bots,prob.terms));
    end
    hangLinks=links(~lands);
    
end