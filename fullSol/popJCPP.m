function [C,C1,C2,C3,C4,C5,CL1,CL2,CL3]=popJCPP(parents,kr,kx,L,nL)

    %efficiency can be improved, if necessary

    a=sqrt(kr./kx);

    iL=(1:nL)';
    sL=nL+iL;
    ssL=nL+iL+1;
    
    nC=parents>0;

    IM=diag(ones(2*nL,1));
    IM(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;
    IM(sub2ind(2*[nL nL],parents(nC),nL+iL(nC)))=-1;

    c0r=-a.*kx.*tanh(a.*L/2);
    c0x=a.*kx./tanh(a.*L);
    c1r=-a.*kx.*tanh(a.*L/2);
    c1x=-a.*kx./sinh(a.*L);
    cSr=2*a.*kx.*tanh(a.*L/2);
    cSx=-a.*kx.*tanh(a.*L/2);

    sFP=[2*nL 2*nL+1];

    FP=zeros(sFP);
    
    %efficiency of sub2ind CAN be improved on with arithmetic operations...
    FP(sub2ind(sFP,iL,parents+1))=c1x;
    FP(sub2ind(sFP,iL,iL+1))=c0x;
    FP(sub2ind(sFP,iL,ssL))=cSx;
    FP(sub2ind(sFP,sL,parents+1))=c1r;
    FP(sub2ind(sFP,sL,iL+1))=c0r;
    FP(sub2ind(sFP,sL,ssL))=cSr;

    C=IM*FP;

    C1=C(iL,1);
    C2=C(iL,1+iL);
    C3=C(iL,ssL);

    CL1=C(sL,1);
    CL2=C(sL,iL+1); %% == CP3'
    CL3=C(sL,ssL); %% == diag(cSr)

    C4=CL2*-inv(C2)*C3+CL3;
    C5=CL2*-inv(C2)*C1+CL1;



end