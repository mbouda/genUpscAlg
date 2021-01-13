function [acroLrs,nxtLr,basiLrs,prvLr]=findDOFs(parents,inLayer,nL)

    tips=setdiff(1:nL,parents);
    acroLrs=cell(nL,1);
    nxtLr=nan(nL,1);
    for t=tips
        l=t;
        p=parents(l);
        layerList=inLayer(l);
        while p>0
            if inLayer(p)~=inLayer(l)
                nxtLr(p)=inLayer(l);
                acroLrs{p}=union(acroLrs{p},layerList);
                layerList=union(layerList,inLayer(p));
            end
            l=p;
            p=parents(l);
        end
    end
    
    %in second loop go down acropetally from base
    bases=find(parents==0)';
    basiLrs=cell(nL,1);
    prvLr=nan(nL,1);
    for b=bases
        l=b;
        %basiLrs{l}=0;  this would be needed if collar boundary condition
        %unknown... but it's known as a condition for solution
        prvLr(l)=0;
        layerList=cat(1,inLayer(l)); %did not add 0 for collar boundary condition
        [basiLrs,prvLr]=basiDOF(l,layerList,basiLrs,prvLr,parents,inLayer); %recursive function to cover all branches
    end
end