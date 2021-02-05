function [eqsOut,nEqs]=remZeroEqs(eqsIn,nEqs)

    cutEq=false(nEqs,1);
    for i=1:nEqs
        cutEq(i)=all(abs(eqsIn(i).coefs)<1e3*eps);
    end

    eqsOut=eqsIn;
    eqsOut(cutEq)=[];
    nEqs=nEqs-sum(cutEq);
end