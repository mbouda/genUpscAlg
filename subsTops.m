function outEqs=subsTops(iLinks,iLayer,inEqs,c1,c2,topLK,inLayer,nTops)

    outEqs=inEqs;
    for i=1:nTops
        eqL=topLK==inLayer(i);
        if inLayer(i)==iLayer
            psiXd=formLinkEq(iLinks(i),sprintf('psiX%d',iLinks(i)),...
                    {sprintf('psi1%d',iLinks(i)),sprintf('G1%d',iLinks(i)),sprintf('psiL%d',inLayer(i))},...
                    [c1(i) c2(i) 1-c1(i)]);
            outEqs(eqL)=subsFor(outEqs(eqL),psiXd.depvar,...
                psiXd.vars,psiXd.coefs);
        else
            varL=ismember(outEqs(eqL).vars,sprintf('psiX%d',iLinks(i)));
            newC=outEqs(eqL).coefs(varL)*cat(2,c1(i),-c2(i),(1-c1(i)));
            newV=cell(1,3);
            newV{1}=sprintf('psi0%d',iLinks(i));
            newV{2}=sprintf('G0%d',iLinks(i));
            newV{3}=sprintf('psiL%d',inLayer(i));

            outEqs(eqL).coefs=cat(2,outEqs(eqL).coefs(~varL),newC);
            outEqs(eqL).vars=cat(2,outEqs(eqL).vars(~varL),newV);
        end
    end

end