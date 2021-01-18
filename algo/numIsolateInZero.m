function  eqOut=numIsolateInZero(eqIn,depvar)


    eqOut=eqIn;
    
    iVar=strcmp(eqIn.vars,depvar);
    coef=-1/eqIn.coefs(iVar);
    
    eqOut.vars(iVar)=[];
    eqOut.coefs(iVar)=[];
    eqOut.depvar=depvar;
    
    eqOut.coefs=coef*eqOut.coefs;
    
    
end