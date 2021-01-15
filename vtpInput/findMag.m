function M=findMag(parents,nL)

    M=zeros(nL,1);
    M(setdiff(1:nL,parents))=1;
    for i=nL:-1:2
        M(parents(i))=M(parents(i))+M(i);
    end

end