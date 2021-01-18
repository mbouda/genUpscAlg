function [elimVars,nVars]=countSysVars(eqs,prob)

    allVars=unique(cat(2,eqs(:).vars));
    domVars=domainBCs(allVars,prob.kLayers);
    elimVars=allVars(~domVars);
    nVars=size(elimVars,2);
end