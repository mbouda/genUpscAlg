function  nDOF=redDOFs(sets,parents,isSub,isX,nSets)

    %counts identical sets
    %for each, says whether there are more of them or they are bigger...
    %minimum of these two values is added to DOF count
        %could probably be improved by accounting for overlaps between
        %non-identical subsets, applying the same logic
    
    %need to prove two crossings are connected on subdomain
    if nSets>1
        
        subN=subNets(parents,isSub,isX,nSets);
        
        catd=false(nSets,1);
        cats=zeros(nSets,1);
        
        iCat=1;
        for i=1:nSets
            if ~catd(i)
                cats(i)=iCat;
                catd(i)=true;
                for j=(i+1):nSets
                    if  subN(i)==subN(j) && ~any(setxor(sets{i},sets{j}))
                        cats(j)=iCat;
                        catd(j)=true;
                    end
                end
                
                %make this into a test for what is less... to re-use later
                iCat=iCat+1;
            end
        end

        nCat=max(cats);
        uqSets=cell(nCat,1);
        lenSet=zeros(nCat,1);
        jCat=zeros(nCat,1);
        oCat=zeros(nCat,1);
        for i=1:nCat
            [~,jCat(i)]=ismember(i,cats);
            oCat(i)=sum(cats==i);
            uqSets{i}=sets{jCat(i)};
            lenSet(i)=size(uqSets{i},1);
        end
        [~,iLen]=sort(lenSet,'ascend');

        soilCat=lenSet<oCat;
        if nCat>1
            for i=1:nCat-1
                iSet=find(iLen==i);
                if soilCat(iSet)
                    j=i+1;  
                    jSet=find(iLen==j);
                    allIn=all(ismember(uqSets{iSet},uqSets{jSet}));
                    sameNet=subN(jCat(iSet))==subN(jCat(jSet));
                    bothSoil=soilCat(jSet);
                    %testing proposition:
                        %all in I are in J
                        %sets are on same subNet
                        %larger category must itself have DOF set by psiS, not nExt 
                    while ~(allIn && sameNet && bothSoil) && j<nCat
                        j=j+1;
                        jSet=find(iLen==j);
                        allIn=all(ismember(uqSets{iSet},uqSets{jSet}));
                        sameNet=subN(jCat(iSet))==subN(jCat(jSet));
                        bothSoil=soilCat(jSet);
                    end

                    if allIn && sameNet && bothSoil
                        %merge i into j
                        lenSet(iSet)=0;   %that is a funcitonal solution, although doesn't really 'merge' fully...
                    end
                end
            end
        end
        nDOF=sum(oCat(~soilCat))+sum(lenSet(soilCat));
    else
        nDOF=nSets;
    end
end