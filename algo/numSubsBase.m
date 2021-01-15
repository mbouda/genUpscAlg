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
            end
            if ismember(sprintf('psi1%d',iDtr),layerEqs(j).vars)
                layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',iDtr))}=sprintf('psi0%d',iSeg);
            end
            if ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars) && ismember(sprintf('G0%d',iSeg),layerEqs(j).vars)
                layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            elseif ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                layerEqs(j)=numPassUpPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
            elseif ismember(sprintf('G0%d',iSeg),layerEqs(j).vars)
                keyboard
                %numPassUpG function is missing for now; unclear if needed
                    %create issue...
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
        dtrLinkI=ismember(prob.iLinks,J);
        if any(termed(dtrLinkI))
            
            j=J(termed(dtrLinkI));
            k=J(~termed(dtrLinkI));
            
            closeEqs(iLinkClose==j).vars{strcmp(closeEqs(iLinkClose==j).vars,...
                    sprintf('psi1%d',j))}=sprintf('psi0%d',iSeg);
            massCons=subsInto(closeEqs(iLinkClose==j),sprintf('G1%d',k),...
                {sprintf('G0%d',iSeg),sprintf('G1%d',j)},...
                [Kx(iSeg)/Kx(k) -Kx(j)/Kx(k)]);

            for j=nLayers:-1:1
                k=J(termed(dtrLinkI));
                if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                    layerEqs(j)=sumVars(layerEqs(j));
                end
                if ismember(sprintf('G1%d',k),layerEqs(j).vars)    
                    layerEqs(j)=subsFor(layerEqs(j),closeEqs(iLinkClose==k).depvar,...
                        closeEqs(iLinkClose==k).vars,...
                        closeEqs(iLinkClose==k).coefs);
                    layerEqs(j)=sumVars(layerEqs(j));
                end
                k=J(~termed(dtrLinkI));
                if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                    layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                end
                if ismember(sprintf('G1%d',k),layerEqs(j).vars)    
                    layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,...
                        massCons.vars,...
                        massCons.coefs);
                    layerEqs(j)=sumVars(layerEqs(j));
                end
                if ismember(sprintf('G0%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                    layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                elseif ismember(sprintf('G0%d',iSeg),layerEqs(j).vars)
                    keyboard  %if ever happens, need to pass up the right stuff here...
                        %create issue
                elseif ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                    layerEqs(j)=numPassUpPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
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
                    for k=J
                        if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                            layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                            layerEqs(j)=sumVars(layerEqs(j));
                        end
                    end
                    if ismember(extraEqs(iExEq).depvar,layerEqs(j).vars)  
                        layerEqs(j)=subsFor(layerEqs(j),extraEqs(iExEq).depvar,extraEqs(iExEq).vars,...
                                                    extraEqs(iExEq).coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(massCons.depvar,layerEqs(j).vars) 
                        layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,massCons.vars,...
                                                    massCons.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end  
                    if ismember(sprintf('G0%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    elseif ismember(sprintf('G0%d',iSeg),layerEqs(j).vars)
                        keyboard  %if ever happens, need to pass up the right stuff here...
                            %create issue
                    elseif ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassUpPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    end
                    if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                        layerEqs(j)=numIsolDep(layerEqs(j));
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
                    for k=J
                        if ismember(sprintf('psi1%d',k),layerEqs(j).vars)
                            layerEqs(j).vars{strcmp(layerEqs(j).vars,sprintf('psi1%d',k))}=sprintf('psi0%d',iSeg);
                            layerEqs(j)=sumVars(layerEqs(j));
                        end
                    end
                    if ismember(massCons.depvar,layerEqs(j).vars) 
                        layerEqs(j)=subsFor(layerEqs(j),massCons.depvar,massCons.vars,...
                                                    massCons.coefs);
                        layerEqs(j)=sumVars(layerEqs(j));
                    end
                    if ismember(sprintf('G0%d',iSeg),layerEqs(j).vars) && ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassUp(layerEqs(j),iSeg,inLayer(iSeg),b2(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    elseif ismember(sprintf('G0%d',iSeg),layerEqs(j).vars)
                        keyboard  %if ever happens, need to pass up the right stuff here...
                            %create issue
                    elseif ismember(sprintf('psi0%d',iSeg),layerEqs(j).vars)
                        layerEqs(j)=numPassUpPsi(layerEqs(j),iSeg,inLayer(iSeg),c1(iSeg),c2(iSeg),c5(iSeg));
                    end
                    if ismember(layerEqs(j).depvar,layerEqs(j).vars)
                        layerEqs(j)=numIsolDep(layerEqs(j));
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