function sol=formSys(eqs,fullSet,iEq)

    notZero=fullSet~=0;
    for i=find(notZero)'
        eqs(i).coefs=cat(2,eqs(i).coefs,-1);
        eqs(i).vars=cat(2,eqs(i).vars,eqs(i).depvar);
        eqs(i).depvar='0';
        eqs(i).kLayer=0;
    end
    
    targEq=eqs(fullSet==iEq);
    remEqs=eqs(fullSet~=iEq);
    nRem=size(remEqs,1); %shows hom many variables we can remove
    nT=nRem+1;
    
    allVars=unique(cat(2,eqs(:).vars));
    N=size(allVars,2);
    
    isPsiX=startsWith(allVars,'psiXBar');
    isPsiL=startsWith(allVars,'psiL');
    isBC=endsWith(allVars,'C');
    isNot=~(isPsiX | isPsiL | isBC);
    
    if sum(isNot)<=nRem
        jL=discoverIndices(allVars(isPsiL),'psiL');
        jX=discoverIndices(allVars(isPsiX),'psiXBar');
        
        j=union(jL,jX);
        dj=j-iEq;
        
        n=sum(isPsiX)+sum(isPsiL);
        domVars=cell(1,n);
        
        
        %prefer:psiX to psiL
        %prefer: lower abs(dj)
        %prefer: -ve over +ve?
        while n>0
            minD=abs(dj)==min(abs(dj));
            if sum(minD)>1
               %need to decide +ve or -ve 
               minD(dj>0)=false;   %remove the positive one, keep negative one
            end
            jMin=j(minD);

            if any(ismember(jX,jMin))
                domVars{n}=sprintf('psiXBar%d',jMin);
                n=n-1;
            end
            if any(ismember(jL,jMin))
                domVars{n}=sprintf('psiL%d',jMin);
                n=n-1;
            end
            
            j(minD)=[];
            dj(minD)=[];
        end
        
        if any(isBC)
            firstL=find(ismember(domVars,sprintf('psiL%d',min(jL))));
            domVars=cat(2,domVars(1:firstL-1),allVars(isBC),domVars(firstL:end));
        end
        varOrd=cat(2,allVars(isNot),domVars);  
        
        cM=zeros(nT,N);
        for i=1:nRem
            for j=1:N
                isVar=ismember(remEqs(i).vars,varOrd{j});
                if any(isVar)
                    cM(i,j)=remEqs(i).coefs(isVar);
                end
            end
        end
        for j=1:N
            isVar=ismember(targEq.vars,varOrd{j});
            if any(isVar)
                cM(nT,j)=targEq.coefs(isVar);
            end
        end
        
        [cM,cutVar]=selectEqs(cM,isNot,nT);
        varOrd(cutVar)=[];
        N=N-sum(cutVar);
        
        nT=size(cM,1);
        nRem=nT-1;
        pres=cM(1:nRem,1:nRem)~=0;
        
        %[cM,pres,nRem,nT]=addRows();
        
        [i,j]=find(pres);
        for k=1:nRem
            I=i(j==k);
            out=cM(k,:);
            if size(I,1)==1
                in=cM(I,:);
                if k~=I
                    cM(k,:)=in;
                    cM(I,:)=out;
                    
                    j(i==I)=[];
                    i(i==I)=[];
                    
                    i(i==k)=I;
                else
                    j(i==I)=[];
                    i(i==I)=[];
                end
                
            elseif size(I,1)>1
                
                ni=size(I,1);
                nI=zeros(ni,1);
                for l=1:ni
                    nI(l)=sum(i==I(l));
                end
                if sum(nI==1)>1
                    keyboard
                    %can't choose... ?
                    %this shows we need to first ID a subset of necessary
                    %equations. Here: the last two lines appear sufficient
                    
                    
                end
                [~,nij]=min(nI);
                I=I(nij);
                
                in=cM(I,:);
                if k~=I
                    cM(k,:)=in;
                    cM(I,:)=out;
                    
                    j(i==I)=[];
                    i(i==I)=[];
                    
                    i(i==k)=I;
                else
                    j(i==I)=[];
                    i(i==I)=[];
                end
            else
        
                if ismember(varOrd{k},eqs(iEq).vars)
                %if present in targEq, then this could be an issue
                 keyboard
                
                else
                    if any(cM(k,nT:N))
                        l=find(cM(k,nT:N),1);
                        inC=cM(:,nRem+l);
                        ouC=cM(:,k);
                        cM(:,k)=inC;
                        cM(:,nRem+l)=ouC;
                    else
                        keyboard
                        %will need to sift through lines k+1:nRem
                        %for ones that have a nonzero entry in column k
                    end
                    
                %else if have another variable that fits this spot 
                    %(preferably psiL, not of iEq)
                    %i.e.: one that has a nonzero coef in row k of cM
                    %then can swith cthose two...
                %can traverse the vars in the order in which they are,
                %since they are already be in order of preference
                        
                    
                %else if was eliminated earlier from eq that is not really
                %needed, can just delete whole Eq...
                
                end
            end
            
        end
        
        for i=1:nRem
            for j=i+1:nT
                cM(j,:)=cM(j,:)-(cM(j,i)/cM(i,i))*cM(i,:);
            end
            cM(i,:)=0;
        end
    else
        keyboard
        %
    end

    mag=abs(cM(end,:));
    hasCoef=mag/max(mag)>1e-16;   %is 16 that too restrictive? set tolerance?
    
    vars=varOrd(hasCoef);
    coefs=cM(end,hasCoef);
    
    depvar=sprintf('psiXBar%d',iEq);
    isDep=ismember(vars,depvar);
    vars(isDep)=[];
    coefs=coefs/(-coefs(isDep));
    coefs(isDep)=[];
    
    sol.depvar=depvar;
    sol.vars=vars;
    sol.coefs=coefs;
    sol.kLayer=iEq;
    
end