function [layerEqs,prob,nLayers]=numConnUDClr(iSeg,layerEqs,lidEqs,extraEqs,iLinkExtra,collarCond,nLayers,prob,termed,Kx,parents)

    par=parents(iSeg);
    
    sibs=parents==par;
    iSibs=find(sibs);
    nSibs=sum(sibs);  %including self
    parLid=find(cat(1,lidEqs.iLink)==par);
    
    if nSibs==1
        switch collarCond
            case 'psiC'
                %lid needs to be passed down to the parent of iSeg 
                    %currently is just still at the collar
                
                lidEqs(parLid).depvar=sprintf('psi1%d',iSeg);
                lidEqs(parLid)=subsFor(lidEqs(parLid),sprintf('G0%d',par),...
                                {sprintf('G1%d',iSeg)},...
                                Kx(iSeg)/Kx(par));
                lidEqs(parLid).iLink=iSeg;
                                
                for j=1:nLayers
                    if ismember(sprintf('psi0%d',par),layerEqs(j).vars)
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),lidEqs(parLid).depvar,...
                                    lidEqs(parLid).vars,...
                                    lidEqs(parLid).coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('G0%d',par),layerEqs(j).vars)
                        layerEqs(j)=subsFor(layerEqs(j),sprintf('G0%d',par),...
                                    {sprintf('G1%d',iSeg)},...
                                    Kx(iSeg)/Kx(par));
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                        layerEqs(j)=numIsolDep(layerEqs(j));
                    end
                end
            case 'G0'
                keyboard
                %need to add substitutions from the lid equations, as per
                %above...
                
                for j=1:nLayers
                    if ismember(sprintf('G0%d',par),layerEqs(j).vars)
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                        layerEqs(j)=subsFor(layerEqs(j),sprintf('G0%d',par),...
                                    {sprintf('G1%d',iSeg)},...
                                    Kx(iSeg)/Kx(par));
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                end
        end
    elseif nSibs==2
        sibLinkI=ismember(prob.iLinks,iSibs);
        sib=setdiff(iSibs,iSeg);
        if any(termed(sibLinkI))
         %if closed: use closeEq.
         %translate this:
             switch collarCond
                case 'psiC'
                keyboard
                    %translate commented analytical to numerical code 
                    %plus need to add propagation of the lid equations...

%                 massCons=isolate(base.botFlux(par),defs.G0(par));
% 
% 
%                 sib=setdiff(find(sibs),iSeg);
%                 base.fclose(sib)=isolate(base.fclose(sib),defs.G1(sib));
% 
% 
%                 for j=prob.kLayers'
%                     base.avg(j)=subsEq(base.avg(j),massCons);
%                     base.avg(j)=subsEq(base.avg(j),base.fclose(sib));
%                     base.avg(j)=subs(base.avg(j),defs.psi1(sib),defs.psi1(iSeg));
%                     base.avg(j)=subs(base.avg(j),defs.psi0(par),defs.psi1(iSeg));
% 
%                 end
                
                case 'GC'
                    keyboard
            end
         
        else
             %if unclosed:
             switch collarCond
                 case 'psiC'
                     if size(cat(1,prob.bots,prob.tops(parents(prob.tops)>0)),1)>nLayers-1
                         %use extraEq
                        iExEq=ismember(iLinkExtra,iSibs);
                        nTargLegs=size(extraEqs(iExEq).targTrack,1);
                        if nTargLegs==2
                            keyboard
                            %how to use each to go both ways?
                        elseif nTargLegs==1
                            keyboard
                            %standard code for using extraEq
                            %also use lidEq

                        end
                     else
                         massCons=formLinkEq(iSeg,sprintf('G0%d',par),...
                            {sprintf('G1%d',iSibs(1)),sprintf('G1%d',iSibs(2))},...
                            [Kx(iSibs(1))/Kx(par) Kx(iSibs(2))/Kx(par)]);
                         massCons=numIsolate(massCons,sprintf('G1%d',sib)); 
                         for j=nLayers:-1:1
                             if ismember(sprintf('G1%d',sib),layerEqs(j).vars) 
                                 %making strong assumption 
                                 layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,massCons.vars,...
                                                            massCons.coefs);
                                 layerEqs(j)=sumVars(layerEqs(j));
                             end
                             if ismember(sprintf('psi1%d',sib),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',sib))}=sprintf('psi1%d',iSeg);
                                layerEqs(j)=sumVars(layerEqs(j));
                             end
                             if ismember(sprintf('psi0%d',par),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi0%d',par))}=sprintf('psi1%d',iSeg);
                                layerEqs(j)=sumVars(layerEqs(j));
                             end 
                             if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                                layerEqs(j)=subsFor(layerEqs(j),sprintf('psi1%d',iSeg),lidEqs(parLid).vars,...
                                                            lidEqs(parLid).coefs);
                                layerEqs(j)=sumVars(layerEqs(j));
                             end
                             if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                                 layerEqs(j)=numIsolDep(layerEqs(j));
                             end
                         end
                     end
                 case  'GC'
                     keyboard
                     %likely should restructure this logic tree before
                     %adding in the GC code...
             end
        end
        
    else
        keyboard
        %for multiple daughters from same junction...
        
    end


end