function [M,mSys]=findMag2(parents,nL)

    M=zeros(nL,1);
    M(setdiff(1:nL,parents))=1;
    for i=nL:-1:2
        if parents(i)>0
            M(parents(i))=M(parents(i))+M(i);
        end
    end
    mSys=sum(M(parents==0));
end
