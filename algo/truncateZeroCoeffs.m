function eqsOut=truncateZeroCoeffs(eqsIn,nEqs)
%     allC=cat(2,eqsIn(:).coefs)';
%     mC=max(abs(allC));

    eqsOut=eqsIn;
%     if any(abs(allC/mC)<1e-16)
        for i=1:nEqs
         %mC=max(abs(eqsIn(i).coefs));     %new line to make this equation-specific
         mC=1; %new line to make this absolute and not relative to other coefs
            del=abs(eqsIn(i).coefs/mC)<1e-16;
            eqsOut(i).vars(del)=[];
            eqsOut(i).coefs(del)=[];
        end
%     end

end