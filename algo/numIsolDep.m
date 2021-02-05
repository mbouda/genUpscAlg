function  eqOut=numIsolDep(eqIn)
%This function isolates the present dependent variable, if it is also
%present on the right hand side
    eqOut=eqIn;
    iVar=strcmp(eqOut.vars,eqIn.depvar);
    if abs(eqIn.coefs(ismember(eqIn.vars,eqIn.depvar))-1)<2*eps
        eqOut.coefs(iVar)=[];
        eqOut.vars(iVar)=[];
        eqOut.depvar='0';
        eqOut.kLayer=0;
    else
        eqOut.coefs=eqOut.coefs(~iVar)/(1-eqOut.coefs(iVar));
        eqOut.vars=eqOut.vars(~iVar);
    end
    
end