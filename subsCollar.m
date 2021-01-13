function  [layerEqs,lidEqs]=subsCollar(iSeg,collarCond,layerEqs,c1,c2,c5,b2,kLayers,inLayer,nSeg)
%for now, making one or other system... depending on bdry condition...

    lidEqs=struct('depvar',cell(nSeg,1),'vars',cell(nSeg,1),'coefs',cell(nSeg,1),'iLink',cell(nSeg,1));
    switch collarCond
        case 'psiC'
            for i=1:nSeg
                psiXf=formLinkEq(iSeg(i),sprintf('psi0%d',iSeg(i)),...
                    {sprintf('psiX%d',iSeg(i)),collarCond,sprintf('psiL%d',inLayer(i))},...
                    [1/c5(i) -1 2-1/c5(i)]);
                psiXb=formLinkEq(iSeg(i),sprintf('psiX%d',iSeg(i)),...
                    {sprintf('psi0%d',iSeg(i)),sprintf('G0%d',iSeg(i)),sprintf('psiL%d',inLayer(i))},...
                    [c1(i) -c2(i) 1-c1(i)]);
                psiXb=subsFor(psiXb,psiXf.depvar,...
                            psiXf.vars,psiXf.coefs);
                psiXb=sumVars(psiXb); 
                psiXb=numIsolDep(psiXb);
                layerEqs(kLayers==inLayer(i))=subsFor(layerEqs(kLayers==inLayer(i)),psiXb.depvar,...
                    psiXb.vars,psiXb.coefs); 
                lidEqs(i)=defLidEqP(iSeg(i),collarCond,c1(i),c2(i),c5(i),inLayer(i)); %gives psi0 as fn of psiC & G0
            end
            
        case 'GC'
            for i=1:nSeg
                psiXa=formLinkEq(iSeg(i),sprintf('G0%d',iSeg(i)),...
                    {sprintf('psiX%d',iSeg(i)),collarCond,sprintf('psiL%d',inLayer(i))},...
                    [-b2(i) 1 b2(i)]);
                psiXb=formLinkEq(iSeg(i),sprintf('psiX%d',iSeg(i)),...
                    {sprintf('psi0%d',iSeg(i)),sprintf('G0%d',iSeg(i)),sprintf('psiL%d',inLayer(i))},...
                    [c1(i) -c2(i) 1-c1(i)]);
                psiXb=subsFor(psiXb,psiXa.depvar,...
                        psiXa.vars,psiXa.coefs);
                psiXb=sumVars(psiXb); 
                psiXb=numIsolDep(psiXb);
                layerEqs(kLayers==inLayer(i))=subsFor(layerEqs(kLayers==inLayer(i)),psiXb.depvar,...
                    psiXb.vars,psiXb.coefs); 
                
                lidEqs=defLidEqsG(iSeg(i),collarCond,c1(i),c2(i),b2(i),inLayer(i)); %gives G0 as fn of GC & psi0
            end
    end
    
    function lidEqs=defLidEqsG(iSeg,collarCond,c1,c2,b2,inLayer)
        keyboard
        %this should be made into a separate file and coded up
        %only here as a placeholder
    end

end