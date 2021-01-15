function targDtr=findTargDtr(pars,parents,prob)

    nPars=size(pars,1);
    targDtr=false(nPars,1);
    for i=1:nPars
        targDtr(i)=any(ismember(find(parents==pars(i)),prob.targ));
    end
end