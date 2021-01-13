function eqOut=subsInto(eqIn,depvar,vars,coefs)
        
    iVar=strcmp(vars,eqIn.depvar);
    eqOut=eqIn;
	eqOut.vars=cat(2,eqIn.vars,vars(~iVar));
    eqOut.coefs=cat(2,coefs(iVar)*eqIn.coefs,coefs(~iVar));
    eqOut.depvar=depvar;            
end