function  eqOut=numIsolDep(eqIn)
%This function isolates the present dependent variable, if it is also
%present on the right hand side

    eqOut=eqIn;
    iVar=strcmp(eqOut.vars,eqIn.depvar);
    eqOut.coefs=eqOut.coefs(~iVar)/(1-eqOut.coefs(iVar));
    eqOut.vars=eqOut.vars(~iVar);

end