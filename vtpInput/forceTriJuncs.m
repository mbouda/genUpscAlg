function newPlant=forceTriJuncs(plant)

    newPlant=plant;
    
    parCount=occur(plant.parents);
    badPar=flipud(unique(plant.parents(parCount>2))); %reversed this in the hope...
        %that doing them in reverse order will prevent having to renumber
        %this array in the loop below...
    
    nBP=size(badPar,1);
    if nBP>0
        warning('Some junctions have >2 daughters, adjusting.','badpar');
    end
    for i=1:nBP

        nD=sum(plant.parents==badPar(i));
        
        if badPar(i)==0 %actually, in this case, want to make maxD ==1?
            maxD=1;  % maximum number of daughters the base is actually allowed
            nNew=nD-maxD; %here, need just one base of network, with parent ==0
            iL=1;
            newU=(0+0.01*(1:nNew))';  %here, assumption that nNew<<100...
        else
            maxD=2;  % maximum number of daughters any segment/link/edge is actually allowed
            nNew=nD-maxD;
            iL=badPar(i);
            newU=fliplr((1-0.01*(1:nNew)))';  %here, assumption that nNew<<100...
        end
        addPlant=breakSegmentAtUs(newU,iL,newPlant,badPar(i));  %in here, currently adds 2 too many dtrs to 1
        
        oldDtrs=find(plant.parents==badPar(i))+nNew;
        dtrs=find(addPlant.parents==badPar(i));
        
            
        badDtrs=dtrs(maxD+1:end);
        newPars=(badPar(i)+(1:nNew))';
        
        for j=1:nNew
            mvPlant=adopt(addPlant,badDtrs(j),newPars(j));
            addPlant=mvPlant;
        end
        newPlant=addPlant;    
    end
        
    newPlant.M=findMag(newPlant.parents,newPlant.nL);

end
