function  [outEqs,iLinkClose]=numCloseTerms(iSeg,b2,c1,c2,jLayer,nTerm)

    CX=b2.*c1./(1-b2.*c2);
    CS=-CX;
    
    CC=[CX CS];
    %want them to be separate outputs by row
    %and then 
    
    %not working vectorised:
    outEqs=struct('coefs',cell(nTerm,1),'vars',cell(nTerm,1),...
        'depvar',cell(nTerm,1));
    iLinkClose=zeros(nTerm,1);
    for i=1:nTerm
        outEqs(i).coefs=CC(i,:);
        outEqs(i).vars={sprintf('psi1%d',iSeg(i)),sprintf('psiL%d',jLayer(i))};
        outEqs(i).depvar=sprintf('G1%d',iSeg(i));
        outEqs(i).iLink=iSeg(i);
        iLinkClose(i)=iSeg(i);
    end
    

end