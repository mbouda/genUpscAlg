function keepExtra=selectExtras(eqs)

    nEqs=size(eqs,1);
    varSet=cell(nEqs,1);
    for i=1:nEqs
        varSet{i}=sort(cat(2,eqs(i).depvar,eqs(i).vars));
    end
    
    cats=zeros(nEqs,1);
    iCat=1;
    for i=1:nEqs
        if cats(i)==0
            cats(i)=iCat;
            for j=i+1:nEqs
                if all(strcmp(varSet{i},varSet{j}))
                    cats(j)=iCat;
                end
            end
            iCat=iCat+1;
        end
    end
    
    keepExtra=false(nEqs,1);
    nCats=iCat-1;
    %now traverse cats
    for iCat=1:nCats
        isCat=cats==iCat;
        nCatEqs=sum(isCat);
        catEqs=eqs(isCat);
        for iEq=1:nCatEqs
            catEqs(iEq)=numIsolate(catEqs(iEq),varSet{find(isCat,1)}(1));
            [catEqs(iEq).vars,iVar]=sort(catEqs(iEq).vars);
            catEqs(iEq).coefs=catEqs(iEq).coefs(iVar);
        end
        
        eqCats=zeros(nCatEqs,1);
        jCat=1;
        for iEq=1:nCatEqs
            if eqCats(iEq)==0
                eqCats(iEq)=jCat;
                for jEq=iEq+1:nCatEqs
                    if all(abs(catEqs(iEq).coefs-catEqs(jEq).coefs)/max(abs(catEqs(iEq).coefs))<1e-12)
                        eqCats(jEq)=jCat;
                    end
                end
                jCat=jCat+1;
            end
        end
        uqEqCats=unique(eqCats);
        keepCat=false(nCatEqs,1);
        for jEqCat=uqEqCats'
            keepCat(find(eqCats==jEqCat,1))=true;
        end
        keepExtra(isCat)=keepCat;
    end
end