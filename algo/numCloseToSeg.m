function [closeEqs,iLinkClose,targTrack]=numCloseToSeg(toSeg,closeEqs,iLinkClose,prob,Kx,b2,c1,c2,c5,parents,inLayer)

    pars=unique(parents(iLinkClose));
    pars=pars(pars~=toSeg);

    %set up list
    targTrack=struct('head',num2cell(iLinkClose),'list',cell(size(iLinkClose)));
    while any(pars)
        
        i=pars(end);

            dtrs=(parents==i);
            J=find(dtrs)';
            if sum(dtrs)==1
                
                clEqI=iLinkClose==J;
                
                
                closeEqs(clEqI).vars{strcmp(closeEqs(clEqI).vars,sprintf('psi1%d',J))}=sprintf('psi0%d',i);
                closeEqs(clEqI).coefs=(Kx(dtrs)/Kx(i))*closeEqs(clEqI).coefs;
                closeEqs(clEqI).depvar=sprintf('G0%d',i);
                closeEqs(clEqI).iLink=i;
                iLinkClose(clEqI)=i;
                          
                closeEqs(clEqI)=numCloseUp(closeEqs(clEqI),i,inLayer(i),c1(i),c2(i),c5(i),b2(i));
           
            elseif sum(dtrs)==2
                %make the bottom flux equation
                massCons=formLinkEq(i,sprintf('G0%d',i),...
                                {sprintf('G1%d',J(1)),sprintf('G1%d',J(2))},...
                                (Kx(J)./Kx(i))');
                for j=J
                    clEqI=iLinkClose==j;
                    massCons=subsFor(massCons,closeEqs(clEqI).depvar,...
                        closeEqs(clEqI).vars,closeEqs(clEqI).coefs);
                    massCons.vars{strcmp(massCons.vars,sprintf('psi1%d',j))}=sprintf('psi0%d',i);
                end
                massCons=sumVars(massCons);
                      
                closeEqs(iLinkClose==J(1))=numCloseUp(massCons,i,inLayer(i),c1(i),c2(i),c5(i),b2(i));
                iLinkClose(iLinkClose==J(1))=i;
                
                closeEqs(iLinkClose==J(2))=[];
                iLinkClose(iLinkClose==J(2))=[];
            else
                keyboard
                %here deal with supernumerary juncitons...
                %although in this case,it's made easier by the fact that we
                %know ALL the daughters are closed
            end
            
            for j=J
                trackHeads=cat(1,targTrack.head);
                isHead=ismember(trackHeads,j);
                if ismember(j,prob.targ)
                    targTrack(isHead).list=cat(1,targTrack(isHead).list,j); %added to list, also pass lists up...
                end
                targTrack(isHead).head=i;
                
            end
            trackHeads=cat(1,targTrack.head);
            isHead=ismember(trackHeads,i);
            if sum(isHead)>1
                newList=cat(1,targTrack(isHead).list);
                targTrack(isHead)=[];
                targTrack(end+1,1).head=i;
                targTrack(end,1).list=newList;
            end
            
        
        pars(end)=parents(i);
        pars=unique(pars);
        pars=pars(pars~=toSeg);
    end
    
    for i=1:size(targTrack,1)
        if ismember(targTrack(i).head,prob.targ)
            targTrack(i).list=cat(1,targTrack(i).list,targTrack(i).head);
        end
    end
    
end