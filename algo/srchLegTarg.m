function [lands,hit]=srchLegTarg(l,parents,targets,bots)

    [lands,j]=ismember(l,targets);
    
    if ~lands
        hit=0;
        if ~ismember(l,bots)
            d=find(parents==l);
            nD=size(d,1);
            i=0;
            while ~lands && i<nD
                i=i+1;
                [lands,hit]=srchTarg(d(i),parents,targets,bots);
            end
        end
    else
        hit=targets(j);
    end


end