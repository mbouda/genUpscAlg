function [closeEqs,layerEqs,termed,iLinkClose]=numCloseInts(prob,closeEqs,iLinkClose,layerEqs,nLayers,Kx,b2,c1,c2,c5,parents,inLayer,termed,jLayer)

    [pars,targDtr]=getPars(termed,prob,parents);
    
    %for each par
    %need to see if any dtr is a prob.targ
    cut=~ismember(inLayer(pars),prob.kLayers) | targDtr; 
    pars(cut)=[];
    
    while any(pars)
  
        i=pars(end);
        
            dtrs=(parents==i);
            if sum(dtrs)==1
                iDtr=find(dtrs);
                clEqI=iLinkClose==iDtr;
                
                %different for each botFlux situation
                closeEqs(clEqI).vars{strcmp(closeEqs(clEqI).vars,sprintf('psi1%d',iDtr))}=sprintf('psi0%d',i);
                closeEqs(clEqI).coefs=(Kx(dtrs)/Kx(i))*closeEqs(clEqI).coefs;
                closeEqs(clEqI).depvar=sprintf('G0%d',i);
                closeEqs(clEqI).iLink=i;
                iLinkClose(clEqI)=i;
                          
                closeEqs(clEqI)=numCloseUp(closeEqs(clEqI),i,inLayer(i),c1(i),c2(i),c5(i),b2(i));
           
                for j=1:nLayers 
                    if ismember(sprintf('G1%d',iDtr),layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),sprintf('G1%d',iDtr),...
                            {sprintf('G0%d',i)},...
                            Kx(i)/Kx(dtrs));
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iDtr))}=sprintf('psi0%d',i);
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=numPassUp(layerEqs(j),i,inLayer(i),b2(i),c1(i),c2(i),c5(i));
                        
%                     elseif ismember(sprintf('psi1%d',iDtr),layerEqs(j).vars)
%                         layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iDtr))}=sprintf('psi0%d',i);
%                         layerEqs(j)=numPassUp(layerEqs(j),i,inLayer(i),b2(i),c1(i),c2(i),c5(i));
                    end
                end

                %placed test here in case the par is actually top, not
                %int... otherwise, unlikely to be needed...
                iLayer=prob.kLayers==inLayer(i);
                if  ismember(sprintf('psiX%d',i),layerEqs(iLayer).vars)
                    layerEqs(iLayer)=subsIntoAvg(layerEqs(iLayer),i,inLayer(i),c1(i),c2(i));
                end
            elseif sum(dtrs)==2

                J=find(dtrs)';

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
                      
                
                iLayer=prob.kLayers==inLayer(i);
                if  ismember(sprintf('psiX%d',i),layerEqs(iLayer).vars)
                    layerEqs(iLayer)=subsIntoAvg(layerEqs(iLayer),i,inLayer(i),c1(i),c2(i));
                end
                
                %unlikely to need to be made more granular as we *know*
                %what is substituted into/from closed segments.
                for j=1:nLayers
                    for k=J
                        if ismember(sprintf('G1%d',k),layerEqs(j).vars)
                            layerEqs(j)=subsFor(layerEqs(j),closeEqs(iLinkClose==k).depvar,...
                                closeEqs(iLinkClose==k).vars,...
                                closeEqs(iLinkClose==k).coefs);
                            layerEqs(j)=sumVars(layerEqs(j));
                            %first sub each closure equation to replace
                            %G1(k) with psi1(k)
                            %then replace psi1(k) with psi0(i)
                            %then sum up coeffs and pass up
                            layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',i);
                            %these last two lines may not be efficient here
                            %but placing outside loop requires separate
                            %test for any(G1) in eq(j)
                            layerEqs(j)=sumVars(layerEqs(j));
                            layerEqs(j)=numPassUpPsi(layerEqs(j),i,inLayer(i),c1(i),c2(i),c5(i));
                        end
                    end
                end
                
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

            %also declare this closed
            termed(ismember(prob.iLinks,i))=true;
        
        if all(termed(ismember(prob.iLinks,find(parents==parents(i))))) && ...
                ~findTargDtr(parents(i),parents,prob) && ...
                ismember(inLayer(parents(i)),prob.kLayers) && ...
                ~ismember(parents(i),pars)
            pars(end)=parents(i); 
        else
            pars(end)=[];
        end
    end
    
end