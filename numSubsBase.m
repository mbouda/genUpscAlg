function [layerEqs,prob,nLayers]=numSubsBase(iSeg,prob,collarCond,closeEqs,iLinkClose,extraEqs,iLinkExtra,layerEqs,nLayers,Kx,b2,c1,c2,c5,termed,parents,inLayer)

% for efficiency, may try put substitution into own average at end, with
% test for necessity or remove entirely ... effect may (or may not) be
% achieved by numPassUp within the loop over all layers

    dtrs=(parents==iSeg);
    if sum(dtrs)==1
          
        iDtr=find(dtrs);
        for j=nLayers:-1:1
            %make Fun for this?
            
            if ismember(sprintf('G1%d',iDtr),layerEqs(j).vars)
                layerEqs(j)=subsFor(layerEqs(j),sprintf('G1%d',iDtr),...
                    {sprintf('G0%d',iSeg)},...
                    Kx(iSeg)/Kx(dtrs));
                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iDtr))}=sprintf('psi0%d',iSeg);
                layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            end
            switch collarCond
                case 'psiC'
                    if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iSeg))}=collarCond;
                    end
                case 'GC'
                    if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                        layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('G1%d',iSeg))}=collarCond;
                    end
            end
            
        end

    elseif sum(dtrs)==2
         
        J=find(dtrs)';
        if any(termed(dtrs))
            j=J(termed(J));
            k=J(~termed(J));
            
            closeEqs(iLinkClose==j).vars{strcmp(closeEqs(iLinkClose==j).vars,...
                    sprintf('psi1%d',j))}=sprintf('psi0%d',iSeg);
            massCons=subsInto(closeEqs(iLinkClose==j),sprintf('G1%d',k),...
                {sprintf('G0%d',iSeg),sprintf('G1%d',j)},...
                [Kx(iSeg)/Kx(k) -Kx(j)/Kx(k)]);

            for j=nLayers:-1:1
                if ismember(sprintf('G1%d',J(1)),layerEqs(j).vars)
                    k=J(termed(J));
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                    layerEqs(j)=subsFor(layerEqs(j),closeEqs(iLinkClose==k).depvar,...
                        closeEqs(iLinkClose==k).vars,...
                        closeEqs(iLinkClose==k).coefs);

                    k=J(~termed(J));
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                    layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,...
                        massCons.vars,...
                        massCons.coefs);

                    layerEqs(j)=sumVars(layerEqs(j));

                    layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                end
                switch collarCond
                    case 'psiC'
                        if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                            layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iSeg))}=collarCond;
                        end
                    case 'GC'
                        if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                            layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('G1%d',iSeg))}=collarCond;
                        end
                end
             end    
        else
            if size(prob.bots,1)>nLayers-1
                iExEq=ismember(iLinkExtra,J);

                massCons=formLinkEq(iSeg,sprintf('G0%d',iSeg),...
                    {sprintf('G1%d',J(1)),sprintf('G1%d',J(2))},...
                    [Kx(J(1))/Kx(iSeg) Kx(J(2))/Kx(iSeg)]);
                massCons=subsFor(massCons,extraEqs(iExEq).depvar,extraEqs(iExEq).vars,...
                    extraEqs(iExEq).coefs);
                massCons=sumVars(massCons);
                massCons=numIsolate(massCons,massCons.vars(startsWith(massCons.vars,'G1')));

                for j=nLayers:-1:1
                    if ismember(sprintf('G1%d',J(1)),layerEqs(j).vars) || ismember(sprintf('G1%d',J(2)),layerEqs(j).vars) 
                        %making strong assumption that J(2) is iLink of extraEq
                        %and both present in layerEq
                        for k=J
                            if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                            end
                        end
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=subsFor(layerEqs(j),extraEqs(iExEq).depvar,extraEqs(iExEq).vars,...
                                                    extraEqs(iExEq).coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,massCons.vars,...
                                                    massCons.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                        if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                            layerEqs(j)=numIsolDep(layerEqs(j));
                        end
                    end
                    switch collarCond
                        case 'psiC'
                            if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iSeg))}=collarCond;
                            end
                        case 'GC'
                            if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('G1%d',iSeg))}=collarCond;
                            end
                    end
                end

                massCons=formLinkEq(iSeg,sprintf('G0%d',iSeg),...
                    {sprintf('G1%d',J(1)),sprintf('G1%d',J(2))},...
                    [Kx(J(1))/Kx(iSeg) Kx(J(2))/Kx(iSeg)]);
                [layerEqs,prob,nLayers]=addLayerEq(extraEqs(iExEq),layerEqs,...
                        prob,nLayers,massCons,iSeg,[],'up',...
                        inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                switch collarCond
                    case 'psiC'
                        if ismember(sprintf('psi1%d',iSeg),layerEqs(end).vars)
                            layerEqs(end).vars{strcmp(layerEqs(end).vars,sprintf('psi1%d',iSeg))}=collarCond;
                        end
                    case 'GC'
                        if ismember(sprintf('G1%d',iSeg),layerEqs(end).vars)
                            layerEqs(end).vars{strcmp(layerEqs(end).vars,sprintf('G1%d',iSeg))}=collarCond;
                        end
                end
            else
                %need to keep one of the G1 variables, as the dof is
                %minimised this way (not sure how general....)
                massCons=formLinkEq(iSeg,sprintf('G0%d',iSeg),...
                    {sprintf('G1%d',J(1)),sprintf('G1%d',J(2))},...
                    [Kx(J(1))/Kx(iSeg) Kx(J(2))/Kx(iSeg)]);
                massCons=numIsolate(massCons,sprintf('G1%d',J(2))); %one is chosen arbitrarily here
                for j=nLayers:-1:1
                    if ismember(sprintf('G1%d',J(2)),layerEqs(j).vars) 
                        %making strong assumption 
                        for k=J
                            if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                            end
                        end
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,massCons.vars,...
                                                    massCons.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                        layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                        if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                            layerEqs(j)=numIsolDep(layerEqs(j));
                        end
                    end
                    switch collarCond
                        case 'psiC'
                            if ismember(sprintf('psi1%d',iSeg),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iSeg))}=collarCond;
                            end
                        case 'GC'
                            if ismember(sprintf('G1%d',iSeg),layerEqs(j).vars)
                                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('G1%d',iSeg))}=collarCond;
                            end
                    end
                end
            end
        end

    else
        keyboard
        %here deal with supernumerary juncitons...
    end
    
end