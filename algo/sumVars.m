function eqOut=sumVars(eqIn)

    % How to do this as efficiently as possible?
    eqOut=eqIn;
    
    uqVars=unique(eqOut.vars);
    nVars=size(uqVars,2);
    for i=1:nVars
        iVar=strcmp(eqOut.vars,uqVars{i});
        if sum(iVar)>1
            eqOut.coefs=cat(2,eqOut.coefs(~iVar),sum(eqOut.coefs(iVar)));
            eqOut.vars=cat(2,eqOut.vars(~iVar),uqVars{i});
        end
    end
    
end