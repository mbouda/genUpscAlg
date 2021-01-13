function eqOut=subsFor(eqIn,var,vars,coefs)

    varI=strcmp(eqIn.vars,var);
    
    eqOut=eqIn;
    eqOut.vars=cat(2,eqIn.vars(~varI),vars);
    eqOut.coefs=cat(2,eqIn.coefs(~varI),eqIn.coefs(varI)*coefs);

end