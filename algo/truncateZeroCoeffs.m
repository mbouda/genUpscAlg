function eqsOut=truncateZeroCoeffs(eqsIn,nEqs)
    
    allC=cat(2,eqsIn(:).coefs)';
    mC=max(allC);
    eqsOut=eqsIn;
    if any(abs(allC/mC)<1e-16)
        for i=1:nEqs
            del=abs(eqsIn(i).coefs/mC)<1e-16;
            eqsOut(i).vars(del)=[];
            eqsOut(i).coefs(del)=[];
        end
    end

end