function U=crossUs3(cz,zLims,inLay,crossLay,nL)

%seems to fail to catch those that go upward...

    U=cell(nL,1);
    for i=1:nL
        if crossLay(i)
            if inLay(i,1)<inLay(i,2)
                l=inLay(i,1):inLay(i,2)-1;
            else
                l=inLay(i,1)-1:-1:inLay(i,2);
            end
            nl=size(l,2);
            for j=1:nl
                U{i}(j)=roots(cz(i,:)-[0 zLims(l(j))]);
            end
        end
        if numel(U{i})>0
            U{i}(U{i}==0 | U{i}==1)=[];
        end
        if size(U{i},2)>1
            U{i}=U{i}';
        end 
    end
    
    
    
    
end