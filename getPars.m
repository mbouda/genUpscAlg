function pars=getPars(termed,prob,parents)

    pars=parents(prob.iLinks(termed)); %parents of termed
    pars=pars(~termed(ismember(prob.iLinks,pars))); %that are open
    nPars=size(pars,1);
    canClose=false(nPars,1);
    for i=1:nPars
        canClose(i)=all(termed(ismember(prob.iLinks,find(parents==pars(i)))));
    end
    pars=unique(pars(canClose));
    
end