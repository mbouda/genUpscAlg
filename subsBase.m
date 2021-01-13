function  layerEqs=subsBase(iSeg,collarCond,layerEqs,c1,c2,kLayers,inLayer,nSeg)
%for now, making one or other system... depending on bdry condition...

%outputs:
%want to substitute psiXd into layerEq
% subs bdry cond into result for appropriate variable



    switch collarCond
        case 'psiC'
            for i=1:nSeg
                psiXd=formLinkEq(iSeg(i),sprintf('psiX%d',iSeg(i)),...
                    {sprintf('psi1%d',iSeg(i)),sprintf('G1%d',iSeg(i)),sprintf('psiL%d',inLayer(i))},...
                    [c1(i) c2(i) 1-c1(i)]);
                psiXd.vars{strcmp(psiXd.vars,sprintf('psi1%d',iSeg(i)))}=collarCond;
                layerEqs(kLayers==inLayer(i))=subsFor(layerEqs(kLayers==inLayer(i)),psiXd.depvar,...
                    psiXd.vars,psiXd.coefs); 
            end
        case 'GC'
            for i=1:nSeg
                psiXd=formLinkEq(iSeg(i),sprintf('psiX%d',iSeg(i)),...
                    {sprintf('psi1%d',iSeg(i)),sprintf('G1%d',iSeg(i)),sprintf('psiL%d',inLayer(i))},...
                    [c1(i) c2(i) 1-c1(i)]);
                psiXd.vars{strcmp(psiXd.vars,sprintf('G1%d',iSeg(i)))}=collarCond;
                layerEqs(kLayers==inLayer(i))=subsFor(layerEqs(kLayers==inLayer(i)),psiXd.depvar,...
                    psiXd.vars,psiXd.coefs); 
            end
    end

end