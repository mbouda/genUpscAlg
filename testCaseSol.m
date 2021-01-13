function fullSol=testCaseSol(testCase)



    sys.parents=testCase.parents;
    sys.inLayer=testCase.inLayer;
    sys.Kr=pi*(2*testCase.params.r-testCase.params.b).*testCase.params.kr./testCase.params.b; 
    sys.Kx=testCase.params.Kx;
    sys.L=testCase.params.L;
    sys.nL=testCase.params.nL;

    psiBC=-0.8e6;
    PSIL=-(0.3:0.2:0.7)*1e6;
    str=strcat('combvec(',repmat('PSIL,',[1 testCase.nDomLayers]),'PSIL)');
    psiBCS=eval(str);
    [nLayer,nBC]=size(psiBCS);
    nC=1;

    fullSol=matSol(sys,psiBC,psiBCS,nBC,nC,nLayer);
    fullSol.psiBCS=psiBCS;
    fullSol.psiBC=psiBC;
    fullSol.nBC=nBC;
    fullSol.nC=nC;

end