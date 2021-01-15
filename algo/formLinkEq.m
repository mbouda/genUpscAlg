function eqOut=formLinkEq(i,depvar,vars,coefs)

    eqOut=struct('depvar',cell(1),'vars',cell(1),'coefs',cell(1),'iLink',cell(1));
    eqOut.depvar=depvar;
    eqOut.vars=vars;
    eqOut.coefs=coefs;
    eqOut.iLink=i;

end