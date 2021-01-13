function [basiLrs,prvLr]=basiDOF(l,layerList,basiLrs,prvLr,parents,inLayer)

        d=parents==l;
        layerList=union(layerList,inLayer(l));
        if any(d)
            crossed=d & inLayer~=inLayer(l);
            if any(crossed)
                for j=find(crossed)'
                    basiLrs{j}=union(basiLrs{j},layerList);
                end
                prvLr(crossed)=inLayer(l);
            end
            
            D=find(d)';
            for j=D
                [basiLrs,prvLr]=basiDOF(j,layerList,basiLrs,prvLr,parents,inLayer);
            end
        end
end