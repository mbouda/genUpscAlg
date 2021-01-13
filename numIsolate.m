function  eqOut=numIsolate(eqIn,depvar)

    eqOut=eqIn;
    
    iVar=strcmp(eqIn.vars,depvar);
    coef=-1/eqIn.coefs(iVar);
    
    eqOut.coefs(iVar)=-1;
    eqOut.vars{iVar}=eqIn.depvar;
    eqOut.coefs=coef*eqOut.coefs;
    eqOut.depvar=depvar;
    
end