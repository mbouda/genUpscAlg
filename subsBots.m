function outEqs=subsBots(iLinks,inEqs,c1,c2,botLK,inLayer,nBots)
    
    outEqs=inEqs;
    for i=1:nBots
        eqL=botLK==inLayer(i);
        varL=ismember(outEqs(eqL).vars,sprintf('psiX%d',iLinks(i)));
        
        newC=outEqs(eqL).coefs(varL)*cat(2,c1(i),c2(i),(1-c1(i)));
        newV=cell(1,3);
        newV{1}=sprintf('psi1%d',iLinks(i));
        newV{2}=sprintf('G1%d',iLinks(i));
        newV{3}=sprintf('psiL%d',inLayer(i));
        
        outEqs(eqL).coefs=cat(2,outEqs(eqL).coefs(~varL),newC);
        outEqs(eqL).vars=cat(2,outEqs(eqL).vars(~varL),newV);
    end

end