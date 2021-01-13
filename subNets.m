function subN=subNets(parents,isSub,isX,nX)

    iX=find(isX)';
    subNets=cell(nX,1);
    onNet=false(nX,1);
    subN=zeros(nX,1);
    subNC=0;
    for i=1:nX
        j=1;
        while ~onNet(i) && j<i
            %check if already member of previous networks...     
            if ismember(iX(i),subNets{j})
                onNet(i)=true;
                subNets{i}=subNets{j};
                subN(i)=subN(j);
            else 
                j=j+1;
            end
        end
        if ~onNet(i)
            subNets{i}=spreadNet(iX(i),iX(i),parents,isSub);
            subNC=subNC+1;
            subN(i)=subNC;
            onNet(i)=true;
        end
    end
    
end